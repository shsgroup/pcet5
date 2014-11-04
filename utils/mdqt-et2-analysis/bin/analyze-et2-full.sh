#!/bin/sh
#
# run analyzis of MDQT trajectories
#

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
      mv $file.tmp $newfile
   done

}


zipfiles() {
   for dist in $*
   do
      gzip $dist
   done
}

analysis_tool="/home/souda/PCET/pcet5/utils/mdqt-et2-analysis/bin/analyze_et2_trajectories_lfs.bin"

curdir=$PWD

for input in $*
do

   inputname=${input%.*}
   dir=$inputname

   dirtraj=`awk -F= '/JOBNAME/{print $2}' ${input}`

   echo -------------------------------------------
   echo Copying trajectories for $inputname ...

   sleep 3

   mkdir -p $curdir/$dir

   if [ -d $curdir/$dir/results ]; then
      echo Removing old data files in $curdir/$dir/results
      rm -rf $curdir/$dir/results/* 2> /dev/null
   else
      mkdir -p $curdir/$dir/results
      echo "creating " $curdir/$dir/results
   fi

   if [ -d $curdir/$dir/trajectories ]; then
      echo Removing old data files in $curdir/$dir/trajectories
      rm -rf $curdir/$dir/trajectories/* 2> /dev/null
   else
      mkdir -p $curdir/$dir/trajectories
      echo "creating " $curdir/$dir/trajectories
   fi

   cp $input $curdir/$dir
   cp $inputname.marcus_parameters $curdir/$dir

   shft=0
   for tarball in `ls ${inputname}_*.tar.*`
   do
      tarext=${tarball##*.}
      if [ $tarext == "bz2" ]; then
         zipoption="j"
      elif [ $tarext == "gz" ]; then
         zipoption="z"
      else
         echo $tarext ": Unknown tarball extension..."
         continue
      fi
      mkdir tmp
      cd tmp
      tar xv${zipoption}f ../$tarball > /dev/null
      cp *.log $curdir/$dir
      cd $dirtraj
      nshift=`ls traj*.dat | awk 'BEGIN{n=0}{n++}END{print n}'`
      echo " ---> " $nshift " trajectories"
      offset $shft `ls traj*.dat`
      cp traj*.dat $curdir/$dir/trajectories
      cp -f ET_free_energy_profiles.dat $curdir/$dir/results
      let "shft+=$nshift"
      cd ../..
      rm -rf tmp
   done
   echo " Total " $shft " trajectory files copied."

   echo
   echo "Now start analysis..."
   echo

   echo "===================================================================================="

   cd $curdir/$dir/results
   echo
   echo "===> Analyzing trajectories in " $curdir/$dir/trajectories
   echo

   trajlist=`ls ../trajectories/trajectory_adiab*.dat 2> /dev/null`
   ntraj=`echo $trajlist | wc -w`

   if [ -s ../trajectories/trajectories.tar.bz2 ] && [ $ntraj == "0" ]; then

      cd ../trajectories
      tar xvjf trajectories.tar.bz2 > /dev/null
      echo "Number of trajectories: " `ls trajectory_adiab*.dat 2> /dev/null | wc -w`
      cd ../results

      eval $analysis_tool ../$inputname.marcus_parameters trajectory_adiab*.dat

      echo -n "===> Archiving binary files with distributions in" $dir/results "."
      zipfiles `ls *_distribution_*.dat`
      echo " Done."
      echo -n "===> Removing files with trajectories in" $dir/trajectories "..."
      rm -rf ../trajectories/trajectory_adiab*.dat
      cd $curdir
      echo " Done."

   elif [ $ntraj -gt 0 ]; then

      echo "Number of trajectories: " $ntraj
      eval $analysis_tool ../$inputname.marcus_parameters ../trajectories/trajectory_adiab*.dat
      echo -n "===> Archiving binary files with distributions in" $dir/results "."
      zipfiles `ls *_distribution_*.dat`
      echo " Done."

      cd $curdir/$dir/trajectories

      if [ ! -e trajectories.tar.bz2 ]; then
         echo -n "===> Archiving files with trajectories in" $dir "..."
         tar cvjf trajectories.tar.bz2 trajectory_adiab*.dat > /dev/null
         echo " Done."
      fi

      echo -n "===> Removing files with trajectories in" $dir "..."
      rm -f trajectory_adiab*.dat
      echo " Done."

      cd $curdir

   else

      echo ":( No trajectories to analyze in " $dir/trajectories

   fi

done

cd $curdir

echo
echo "========== ALL DONE ================================================================"
echo
