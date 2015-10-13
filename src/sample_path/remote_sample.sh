#!bin/bash

echo Connecting to the server
ssh giverny


echo Initiate R script
Rscript --vanilla sample_path_script.R &

disown
logout

EOF

