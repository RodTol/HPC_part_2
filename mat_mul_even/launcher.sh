#!/bin/bash
echo "Welcome to the launcher for the matrix multiplication"
echo "This program works only for evenly distributed matrix!"
read -p "Insert the matrix dimension N: " N
read -p "Insert the number of processor: " n_proc

while true; do
    read -p "Do you wish to run in debug mode? " yn
    case $yn in
        [Yy]* ) make run flags=debug N=$N CORES=$n_proc;break;;
        [Nn]* ) make run N=$N CORES=$n_proc;break;;
        * ) echo "Please answer yes or no.";;
    esac

done
