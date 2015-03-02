#!bin/bash

ssh giverny "
  cd ~/Masters/src/predictions/behavior
  nohup Rscript --vanilla qp_plot_script.R > /dev/null 2>&1 &
"
