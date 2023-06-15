#!/bin/bash

salloc --exclusive -A tra23_units_bis -p m100_usr_prod --time=2:0:0 \
 -N 1 \
--ntasks-per-node=32 \
 --gres=gpu:4 --mem=246000 --ntasks-per-core=1 \
--job-name=tolloi_mat_mul
