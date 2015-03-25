
#!bin/bash

ssh giverny "
  cd ~/Masters/src/fitting/
  nohup Rscript --vanilla validate_sampler_script.R > /dev/null 2>&1 &
"