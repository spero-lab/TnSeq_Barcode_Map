#!/bin/bash

# Changing file permissions for all .sh scripts in scripts folder (in case not already done) 
chmod 755 ./* 

# Running Scripts 
    # Set up Conda Env  
    ./setup1_conda.sh
    # Cp Input Files 
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate TNSeq_env
    ./setup2_organize.sh
exit