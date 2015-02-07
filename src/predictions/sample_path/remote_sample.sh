#!bin/bash

ssh giverny "
  cd ~/Masters/src/predictions/sample_path
  nohup Rscript --vanilla sample_path_script.R > /dev/null 2>&1 &
"
