#include <stdio.h>
#include <stdlib.h>
#include "headers/utilities.h"

int FileExists(const char *filename)
{    
   FILE *fp = fopen (filename, "r");
   if (fp!=NULL) { fclose (fp); };
   return (fp!=NULL);
}


/**
 * @brief This function print on a file the data. 
 * NOTE: this function is designed for data distributed on the n1 direction
 * 
 * @param name name of the file
 * @param n1 global grid dimension 
 * @param n2 global grid dimension
 * @param n3 global grid dimension
 * @param n1_local local grid dimension for the n1 dimension
 * @param n1_local_offset local offeset for the n1 dimension
 * @param dir direction of the print
 * @param data the data array
 */
void plot_data_2d( char* name, int n1, int n2, int n3, int n1_local, int  n1_local_offset, int dir, double* data )
{
    int i1, i2, i3, i;
    FILE *fp;
    int num = 1;    
    char buf[256];
    int index;

    int mype, npes, owner;
    int *sizes, *displ;
    double* buffer, *buffer1d, *local_buffer;

    MPI_Comm_rank( MPI_COMM_WORLD, &mype );
    MPI_Comm_size( MPI_COMM_WORLD, &npes );

    snprintf(buf, sizeof(buf), "%s.dat", name); 

    owner=npes+1;
    if ( (n1/2 > n1_local_offset) &&  (n1/2 <= n1_local_offset + n1_local) ) {
        owner = mype;
    }

    /*
    * HINT: Assuming you sliced your system along i3, iif idir==1 or idir==2 you 
    *       need to take the correct slice of the plane from each process. If idir==3, you need  
    *       to understand which process holds the plane you are inrested in. 
    */
    if ( dir == 1)
        {
        i1=n1/2-1; 
        if ( mype == owner)
            {
            fp = fopen (buf, "w");
            for (i2 = 0; i2 < n2; ++i2)
                {
                for (i3 = 0; i3 < n3; ++i3)
                    {
                    index = index_f(i1-n1_local_offset,i2,i3,n1_local,n2,n3);  
                    fprintf(fp, " %14.6f ", data[index] );
                    }
                fprintf(fp, "\n");
                }
            fclose(fp);
            }
        }  
    else if ( dir == 2) {
        i2=n2/2-1;
        sizes = (int*)malloc(npes*sizeof(int));
        displ = (int*)calloc(npes,sizeof(int));

        buffer= (double*)malloc(n1*n3*sizeof(double));
        buffer1d= (double*)malloc(n1*sizeof(double));
        local_buffer = (double*)malloc(n1_local*sizeof(double));

        MPI_Gather(&n1_local, 1, MPI_INT, sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if ( mype  ==  0) { /*Solo il master ha la matrice completa*/
            for (i=1; i < npes; ++i) {
                displ[i] = sizes[i-1] + displ[i-1];
            }
        }

        for (i3=0; i3< n3; ++i3) {
            for ( i1 = 0; i1 < n1_local; ++i1) {
                index = index_f (i1,i2,i3, n1_local, n2 , n3);
                local_buffer[i1] = data[index];   
            }                
            MPI_Gatherv(local_buffer,n1_local, MPI_DOUBLE,buffer1d,sizes, displ, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            for ( i1 = 0; i1 < n1; ++i1) {
                buffer[ i1*n3 + i3] = buffer1d[i1];
            }
        }

        if (mype == 0) { 
            fp = fopen (buf, "w");
            for (i1 = 0; i1 < n1; ++i1) {
                for (i3 = 0; i3 < n3; ++i3) {
                    fprintf(fp, " %14.6f ", buffer[i1*n3 + i3] );
                }
                fprintf(fp, "\n");
            }
            fclose(fp);
        }

        free(sizes);
        free(displ);
        free(buffer);
        free(buffer1d);
        free(local_buffer);
        
    } else if ( dir == 3) {
        i3=n3/2-1;
        sizes = (int*)malloc(npes*sizeof(int));
        displ = (int*)calloc(npes,sizeof(int));

        buffer= (double*)malloc(n1*n2*sizeof(double));
        buffer1d= (double*)malloc(n1*sizeof(double));
        local_buffer = (double*)malloc(n1_local*sizeof(double));

        MPI_Gather(&n1_local, 1, MPI_INT, sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if ( mype  ==  0) { 
            for (i=1; i < npes; ++i) {
                displ[i] = sizes[i-1] + displ[i-1];
            }        
        }
        for (i2=0; i2 < n2; ++i2) {
            for ( i1 = 0; i1 < n1_local; ++i1) {
                index = index_f (i1,i2,i3, n1_local, n2 , n3);
                local_buffer[i1] = data[index];
            }
            MPI_Gatherv(local_buffer,n1_local, MPI_DOUBLE,buffer1d,sizes, displ, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            for ( i1 = 0; i1 < n1; ++i1) {
                buffer[ i1*n2 + i2] = buffer1d[i1];
            }
        }
        if (mype == 0) { 
            fp = fopen (buf, "w");
            for (i1 = 0; i1 < n1; ++i1) {
                for (i2 = 0; i2 < n2; ++i2) {
                    fprintf(fp, " %14.6f ", buffer[i1*n2 + i2] );
                }
                fprintf(fp, "\n");
            }
            fclose(fp);
        }

        free(sizes);
        free(displ);
        free(buffer);
        free(buffer1d);
        free(local_buffer);
    } else {
        fprintf(stderr, " Wrong value for argument 7 in plot_data_2d \n");
    }
}
