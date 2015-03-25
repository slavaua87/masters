#!bin/bash

ssh giverny "
  cd ~/Masters/src/fitting/
  nohup Rscript --vanilla posterior_sampler_script.R > /dev/null 2>&1 &
"
