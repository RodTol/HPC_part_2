<img src="images/cineca.png" 
        alt="Picture" 
        width="400" 
        height="100" 
        style="display: block; margin: 0 auto" />
# Advanced parallel programming
#### Rodolfo Tolloi
### Built with: 
![<C>](https://img.shields.io/badge/C-00599C?style=for-the-badge&logo=c&logoColor=white)
![nVIDIA](https://img.shields.io/badge/nVIDIA-%2376B900.svg?style=for-the-badge&logo=nVIDIA&logoColor=white)  
Homeworks and exercises from the course in Advanced Parallel Programming.
 Each exercise has its own directory and inside of it you'll find the code and the results of the benchmarking (if requested).  
### TODO list:
1. Print function for all matrices duplicates data, so the max size for a matrix in a node is half the memory of the node  
~~2. My initialisation of the matrices is done by working on the local buffer of each matrix. This can be a problem if we want to specify a generic rule for the whole A matrix. For example if I want the identity. I need to create a rule using the global indexes and the offsets (I done this for the identity and jacobi).~~
3. Fix comments for all the exercises