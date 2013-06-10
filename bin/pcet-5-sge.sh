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
SGEQUEUE="sing*.q"
WALLTIME="720:00:00"
COMP="intel"
SCRBASE="local"
THREADS=1
NJOBS=1
inputlist=""

#-- help
if [ x$1 == x ]; then
    echo "================================================================================="
    echo "Synopsis:"
    echo $0 "<options> <input_file_name>"
    echo "---------------------------------------------------------------------------------"
    echo "Options:"
    echo "-q | --queue         :  SGE queue (default batch)"
    echo "-w | --walltime      :  wall time (default 240 hours)"
    echo "-s | --scratch-base  :  global/local - base scratch directory (default local)"
    echo "-c | --compiler      :  intel/pgi - build compiler (default intel)"
    echo "-t | --threads       :  number of OpenMP threads used by diagonalization routines"
    echo "-j | --jobs          :  number of jobs to submit (default 1)"
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
   /usr/bin/modulecmd bash load intel/intel-11
   if [ "$THREADS" -gt "1" ]; then
      EXEC=$BINDIR/pcet_5.3_intel_11.1.omp.x
   else
      EXEC=$BINDIR/pcet_5.3_intel_11.1.x
   fi
   /usr/bin/modulecmd bash list
elif [ x$COMP == xpgi ]; then
   /usr/bin/modulecmd bash load pgi64
   /usr/bin/modulecmd bash list
   EXEC=$BINDIR/pcet_5.3_pgi-9.02-x86_64.x
else
   echo "Unknown build version: " $COMP
   exit 1
fi


if [ x$SCRBASE == xglobal ]; then
   echo Scratch directory - global scratch space on the master node: /scratch/${USER}
   SCRATCH=/scratch
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
module load openmpi/intel-11

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

# deliver the output directory
tar cvjf $JOBNAME.tar.bz2 *
cp $JOBNAME.tar.bz2 \$SGE_O_WORKDIR

echo "Removing SGE scripts..."
rm -f \$SGE_O_WORKDIR/*.sge

echo "Removing *po files..."
rm -f \$SGE_O_WORKDIR/*.po*

echo "Removing scratch directories..."
rm -rf \$SCRDIR

module purge

EnD

qsub ${HERE}/${JOBNAME}.sge
#sleep 1.2

done

echo
echo `date` : ${NJOBS} jobs has been submitted to ${SGEQUEUE} queue

done

exit 0
