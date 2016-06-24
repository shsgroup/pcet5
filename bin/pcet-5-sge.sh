#!/bin/sh
#
# submits PCET'5 job to the SGE queue on shscluster@SCS.UIUC
# Synopsis:
# <name-of-the-script> <options>
#
# Alexander Soudackov, University of Illinois at Urbana-Champaign, 02/13/2013
#

HERE=`pwd`

#-- directory where pcet executables reside (change it if needed)
BINDIR=/home/souda/PCET/pcet5/bin

#-- defaults
SGEQUEUE="single.q"
WALLTIME="168:00:00"
COMP="intel"
SCRBASE="local"
THREADS=1
NJOBS=1
EMAIL=""
VERSION="5.4"
inputlist=""

#-- help
if [ x$1 == x ]; then
    echo "================================================================================="
    echo "Synopsis:"
    echo $0 "<options> <input_file_name>"
    echo "---------------------------------------------------------------------------------"
    echo "Options:"
    echo "-q | --queue         :  SGE queue (default single.q)"
    echo "-w | --walltime      :  wall time - default 168 hours (7 days)"
    echo "-s | --scratch-base  :  global/global2/local - base scratch directory (default local)"
    echo "-c | --compiler      :  intel/pgi - build compiler (default intel)"
    echo "-t | --threads       :  number of OpenMP threads used by diagonalization routines"
    echo "-j | --jobs          :  number of jobs to submit (default 1)"
    echo "-m | --mail          :  send mail upon successful job completion"
    echo "-v | --version       :  program version (default 5.4)"
    echo "================================================================================="
    exit 0
fi

#-- read options
while test "x$1" != x ; do
case $1 in
-q | --queue        ) shift; SGEQUEUE="$1" ;;
-w | --walltime     ) shift; WALLTIME="$1" ;;
-s | --scratch-base ) shift; SCRBASE="$1" ;;
-c | --compiler     ) shift; COMP="$1" ;;
-t | --threads      ) shift; THREADS="$1" ;;
-j | --jobs         ) shift; NJOBS="$1" ;;
-m | --mail         ) shift; EMAIL="$1" ;;
-v | --version      ) shift; VERSION="$1" ;;
*) break ;;
esac
shift
done

inputlist=$*

if [ "x$inputlist" == "x" ]; then
    echo "you MUST specify at least one <input_file_name> "
    exit 1
fi


if [ x$COMP == xintel ]; then
   module load intel/intel-11
   if [ "$THREADS" -gt "1" ]; then
      EXEC=$BINDIR/pcet_${VERSION}_intel_11.1.omp.x
   else
      EXEC=$BINDIR/pcet_${VERSION}_intel_11.1.x
   fi
   module list
elif [ x$COMP == xpgi ]; then
   module load pgi
   module list
   EXEC=$BINDIR/pcet_${VERSION}_pgi-13.9-x86_64.x
elif [ x$COMP == xopen64 ]; then
   module load open64
   module load acml/open64
   module list
   EXEC=$BINDIR/pcet_${VERSION}_open64-x86_64.x
else
   echo "Unknown build version: " $COMP
   exit 1
fi


if [ x$SCRBASE == xglobal ]; then
   echo Scratch directory - global scratch space on the master node: /share/scratch/${USER}
   SCRATCH=/share/scratch
elif [ x$SCRBASE == xglobal2 ]; then
   echo Scratch directory - global scratch2 space on the master node: /share/scratch2/${USER}
   SCRATCH=/share/scratch
elif [ x$SCRBASE == xlocal ]; then
   echo Scratch directory - local scratch space on the compute node: /scr/${USER}
   SCRATCH=/state/partition1
else
   echo Default scratch directory - local scratch space on the compute node: /scr/${USER}
   SCRATCH=/state/partition1
fi


for INPUTFILE in $inputlist
do

if [ ! -e $INPUTFILE ]; then
    echo "Input file " $INPUTFILE " does not exist :("
    break
fi

INPUTNAME=${INPUTFILE%.*}
LOGFILE=${INPUTNAME}.log

echo
echo "----------------------------------------"
echo "Input parameters:"
echo "----------------------------------------"
echo "INPUTFILE = " $INPUTFILE
echo "LOGFILE   = " $LOGFILE
echo "INPUTNAME = " $INPUTNAME
echo "THREADS   = " $THREADS
echo "SCRTYPE   = " $SCRBASE
echo "SGEQUEUE  = " $SGEQUEUE
echo "NJOBS     = " $NJOBS
if [ x$EMAIL != x ]; then
   echo "----------------------------------------"
   echo "Notification will be sent to " $EMAIL
fi
echo "----------------------------------------"
echo


for job in `seq 1 $NJOBS`
do

JOBNAME=${INPUTNAME}_$$_job_$job
JOBDIR=${INPUTNAME}_$$_job_$job

#########################
# create the SGE script
#########################

cat << EnD > $JOBNAME.sge
#$ -S /bin/sh
#$ -N $JOBNAME
#$ -q $SGEQUEUE
#$ -l h_rt=$WALLTIME
#$ -j yes
#$ -cwd
#$ -V
#$ -pe mpi1 $THREADS

module load intel/intel-11
module load openmpi/1.6.5/intel-11
module load pgi
module load open64
module load acml/open64

# define and create the scratch directories
SCRDIR=$SCRATCH/\${SGE_O_LOGNAME}/\${JOB_NAME}_\${JOB_ID}
mkdir \$SCRDIR

cd \$SGE_O_WORKDIR

# extract the name of output directory and the default name of the checkpoint file
OUTPUTDIR=`awk -v IGNORECASE=1 '/JOBNAME/{print substr(\$0,index(\$0,"=")+1)}' $INPUTFILE`

# extract the names of additional input files

ADDINP=\`awk -v IGNORECASE=1 '/GEOM/  {i1=index(\$0,"GEOM=" )+5;i2=index(substr(\$0,i1),")");split(substr(\$0,i1,i2-1),a,",");print a[1]}
                              /PARS/  {i1=index(\$0,"PARS=" )+5;i2=index(substr(\$0,i1),")");split(substr(\$0,i1,i2-1),a,",");print a[1]}
                              /TREAD/ {i1=index(\$0,"TREAD=")+6;i2=index(substr(\$0,i1),")");split(substr(\$0,i1,i2-1),a,",");print a[1]}
                              /DKLIN/ {i1=index(\$0,"DKLIN=")+6;i2=index(substr(\$0,i1),")");split(substr(\$0,i1,i2-1),a,",");print a[1]}' $INPUTFILE\`

DCHKFILE=\`awk -v IGNORECASE=1 '/RESTART=/{i1=index(\$0,"RESTART=" )+8;i2=index(substr(\$0,i1),")");split(substr(\$0,i1,i2-1),a,",");print a[1]}' $INPUTFILE\`

if [ x\$DCHKFILE == x ]; then
   LINE=\`grep RESTART $INPUTFILE\`
   if [ x\$LINE != x ]; then
      DCHKFILE=\$OUTPUTDIR.dchk
   fi
fi

ADDINP="\$ADDINP \$DCHKFILE"

echo "Additional input files: " \$ADDINP

# copy input file and checkpoint file to the scratch (execution) directory
if [ -e $INPUTFILE ] ; then
    cp ${INPUTFILE} \$SCRDIR
fi

for file in \$ADDINP
do
   if [ -e \$file ]; then
        cp -f \$file \$SCRDIR
   else
      echo "File " \$file " is referenced in the main input but can not be found."
      echo "Check your input..."
      echo "Removing scratch directories..."
      rm -rf \$SCRDIR
      exit 1
   fi
done

cd \$SCRDIR

#NODELIST=\`cat \$PBS_NODEFILE | tr -cs "[:alnum:]" "," | rev | cut -b 2- | rev\`
echo "=========================================================================="
echo "Starting PCET on " \$NSLOTS " SGE slots"
echo "Date: " \`date\`
echo "--------------------------------------------------------------------------"
echo "Number of OpenMP threads: " $THREADS
echo "--------------------------------------------------------------------------"
echo "SCRATCH directory: " \$SCRDIR
echo "=========================================================================="

export OMP_NUM_THREADS=$THREADS
export MKL_NUM_THREADS=$THREADS
time $EXEC $INPUTFILE | tee $LOGFILE

echo "==========================================================================="
echo "Job is finished."
echo "---------------------------------------------------------------------------"
echo "Date: " \`date\`
echo "==========================================================================="

# send e-mail about the job completion
if [ x$EMAIL != x ]; then
   MAILTEXT=msg$RANDOM
   rm -rf \$MAILTEXT
   echo "Subject: The job $JOBNAME has finished"  >> \$MAILTEXT
   echo "To: $EMAIL"  >> \$MAILTEXT
   echo "From: souda@$HOSTNAME"  >> \$MAILTEXT
   echo "Reply-To: souda@$HOSTNAME"  >> \$MAILTEXT
   echo "$HOSTNAME : The job $JOBNAME has finished" >> \$MAILTEXT
   cat \$MAILTEXT | sendmail -t
   rm \$MAILTEXT
fi

# deliver the output directory
tar cvjf $JOBNAME.tar.bz2 *
cp $JOBNAME.tar.bz2 \$SGE_O_WORKDIR

echo "Removing SGE scripts..."
rm -f \$SGE_O_WORKDIR/*.sge

echo "Removing *po files..."
rm -f \$SGE_O_WORKDIR/*.po*

echo "Removing scratch directories..."
rm -rf \$SCRDIR

EnD

qsub ${HERE}/${JOBNAME}.sge
#sleep 1.2

done

echo
echo `date` : ${NJOBS} jobs has been submitted to ${SGEQUEUE} queue

done

exit 0
