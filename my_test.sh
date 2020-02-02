#!/bin/bash

docker run -it -v $pwd:/data fegenie ./FeGenie.py -bin_dir my_test_dataset -bin_ext txt -out fegenie_out -hmm_lib HMM-lib -t $(nproc)
