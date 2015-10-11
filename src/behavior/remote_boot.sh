
#!bin/bash

ssh giverny "
  cd ~/Masters/src/predictions/behavior
  nohup Rscript --vanilla test_boot.R > /dev/null 2>&1 &
"
