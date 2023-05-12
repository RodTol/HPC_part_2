#!/bin/bash

salloc -A tra23_units -p m100_usr_prod --time=1:0:0 \
 -N 2 \
--ntasks-per-node=32 \
 --gres=gpu:4 --mem=246000 --ntasks-per-core=1 \
--job-name=tolloi_mat_mul 
