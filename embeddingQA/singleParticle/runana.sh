#! /bin/bash

array=( pip pim kp km ap p )

for i in "${array[@]}"; do
  particle=$i
  echo processing $particle
  intputfile="${particle}/production_${particle}.root"
  outputfile="out_${particle}.root"
  root -b -q readtree.C\(\"${intputfile}\",\"${outputfile}\"\)
done

