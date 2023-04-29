#!/bin/bash
echo "Welcome to the launcher for the matrix multiplication"
echo "This program works only for evenly distributed matrix!"
read -p "Insert the matrix dimension N: " N
read -p "Insert the number of processor: " n_proc

while true; do
    read -p "Do you wish to run in debug mode? " yn
    case $yn in
        [Yy]* ) FLAGS='debug';break;;
        [Nn]* ) echo "NO DEBUG";break;;
        * ) echo "Please answer yes or no.";;
    esac

done

while true; do
    echo "What kind of compilation you want?"
    echo "- Standard MPI (0)"
    echo "- OpenBlas DGEMM  (1)"
    echo "- cuBLAS DGEMM (2)"
    read -p "select one:" input
    case $input in
        [0]* ) COMPILATION='';break;;
        [1]* ) COMPILATION='dgemm';break;;
        [2]* ) COMPILATION='gpu';break;;
        * ) echo "Please answer 1 2 or 3.";;
    esac

done

make run flags=$FLAGS compilation=$COMPILATION N=$N CORES=$n_proc
