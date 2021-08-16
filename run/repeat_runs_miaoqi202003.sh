#!/bin/bash

for jj in $(seq 1 1 20);	
do {	


  # echo "miaoqi202003_G0029_star_run; quit" | matlab -nodisplay
  echo "miaoqi202003_G0029_star_run; quit" | ~/Documents/MATLABR2019B/bin/matlab -nodisplay

  # echo "miaoqi202003_G0024_apslogo_run; quit" | matlab -nodisplay 
  # echo "miaoqi202003_G0024_apslogo_run; quit" | ~/Documents/MATLABR2019B/bin/matlab -nodisplay


  var="$jj""-run"
  mkdir $var
    	
  mv "G0029_star.mat" ./$var
  # mv "G0024_apslogo.mat" ./$var

  mv *.jpg ./$var
} done



