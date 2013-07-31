#!/bin/sh

offset() {

   off=$1
   shift
   files=$*

   for file in $files
   do
      mv $file $file.tmp
   done

   for file in $files
   do
      num=${file:17:6}
      newnum=`echo $num "+" $off + 1000000|bc`
      if [ "$newnum" -ge "0" ]; then
         newnum=${newnum:1}
      else
	 newnum=$num
      fi
      newfile=${file/$num/$newnum}
      #echo -n $file " ---> " $newfile
      mv $file.tmp $newfile
      #echo "   ok."
   done

}

#scrdir=/scratch/souda/mdqt-et2-analysis/et-o2c-er27.6251
scrdir=$PWD
#shift
mkdir -p $scrdir

for input in $*
do

   inputname=${input%.*}
   dirtraj=`awk -F= '/JOBNAME/{print $2}' ${input}`
   echo -------------------------------------------
   echo Copying $inputname ...
   mkdir $scrdir/$inputname
   mkdir $scrdir/$inputname/results
   cp $input $scrdir/$inputname

   shft=0
   for tarball in `ls ${inputname}_*.tar.bz2`
   do
      mkdir tmp
      cd tmp
      tar xvjf ../$tarball > /dev/null
      cp *.log $scrdir/$inputname
      cd $dirtraj
      nshift=`ls traj*.dat | awk 'BEGIN{n=0}{n++}END{print n}'`
      echo " ---> " $nshift " trajectories"
      offset $shft `ls traj*.dat`
      cp traj*.dat $scrdir/$inputname
      let "shft+=$nshift"
      cd ../..
      rm -rf tmp
   done
   echo " Total " $shft " trajectory files copied."

done

echo
echo "All done!"
echo

