#!/bin/bash

input=$1
output=$2

acts-install/bin/ActsAnalysisMaterialComposition \
      -i $input \
      -o $output \
      --sub-names beampipe pixel sstrips lstrips solenoid \
      --sub-rmin 0:25:202:720:1160 \
      --sub-rmax 24.4:200:720:1140:1200 \
      --sub-zmin -4000:-2400:-3150:-3150:-3000 \
      --sub-zmax 4000:2400:3150:3150:3000 \
      -s
