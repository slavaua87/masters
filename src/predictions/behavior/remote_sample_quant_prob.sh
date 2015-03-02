#!bin/bash

ssh giverny "
  cd ~/Masters/src/predictions/behavior
  nohup Rscript --vanilla quant_prob_script.R > /dev/null 2>&1 &
"
