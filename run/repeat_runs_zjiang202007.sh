#!/bin/bash

export LD_LIBRARY_PATH=~/Documents/MATLAB/HDF5Plugin
export HDF5_PLUGIN_PATH=~/Documents/MATLAB/HDF5Plugin


for jj in $(seq 1 1 10);	
do {	


  # echo "zjiang202007_L0367_A4rot2Dptycho000deg_run; quit" | matlab -nodisplay

  echo "zjiang202010_ptycho0080star_run; quit" | ~/Documents/MATLABR2019B/bin/matlab -nodisplay
  # echo "zjiang202010_ptycho0078apslogo_run; quit" | ~/Documents/MATLABR2019B/bin/matlab -nodisplay
  # echo "zjiang202007_L0159_petraIII_part0_run; quit" | ~/Documents/MATLABR2019B/bin/matlab -nodisplay
  # echo "zjiang202007_L0160_petraIII_part1_run; quit" | ~/Documents/MATLABR2019B/bin/matlab -nodisplay
  # echo "zjiang202007_L0161_petraIII_part2_run; quit" | ~/Documents/MATLABR2019B/bin/matlab -nodisplay

  var="$jj""-run"
  mkdir $var



  mv "ptycho0080_star_125nm.mat" ./$var  
  # mv "ptycho0078_apslogo_125nm.mat" ./$var  
  # mv "L0159_petraIII_part0.mat" ./$var
  # mv "L0160_petraIII_part1.mat" ./$var
  # mv "L0161_petraIII_part2.mat" ./$var

  mv *.jpg ./$var

} done



