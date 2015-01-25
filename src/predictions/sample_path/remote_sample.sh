#!bin/bash

echo Connecting
ssh vyacheslav@giverny <<EOF

pwd

echo Initiate R script
Rscript sample_path_script.R

EOF

