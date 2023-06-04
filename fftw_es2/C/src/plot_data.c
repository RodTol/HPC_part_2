/*
 * Assignment 
 * 
 * Adapt these routines to handle distributed arrays 
 * 
 * Created by G.P. Brandino, I. Girotto, R. Gebauer
 * Last revision: March 2016
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "headers/utilities.h"

int FileExists(const char *filename)
{    
   FILE *fp = fopen (filename, "r");
   if (fp!=NULL) fclose (fp);
   return (fp!=NULL);
}


void plot_data_1d( char* name, int n1, int n2, int n3, int dir, double* data)
{
    int i1, i2, i3;
    FILE *fp;
    int num = 1;    
    char buf[256];
    int index;

    snprintf(buf, sizeof(buf), "%s_%d.dat", name, num); 
    while (FileExists(buf))
          {
          num++;
          snprintf(buf, sizeof(buf), "%s_%d.dat", name, num);
          }
    fp = fopen (buf, "w");

 /*
  * HINT: Assuming you sliced your system along i3, iif idir.eq.1 or idir.eq.2 you 
  *       need to understand which process holds the line you are inrested in. If idir.eq.3
  *       you need to take the correct slice of the line from each process
  */

    if ( dir == 1)
        {
        i2=n2/2-1;
        i3=n3/2-1;
        for (i1 = 0; i1 < n1; ++i1)
            {
            index = index_f(i1,i2,i3,n1,n2,n3);  
            fprintf(fp, " %14.6f ", data[index] );
            }
        fprintf(fp, "\n");
        }
    else if ( dir == 2)
        {
        i1=n1/2-1;
        i3=n3/2-1;
        for (i2 = 0; i2 < n2; ++i2)
            {
            index = index_f(i1,i2,i3,n1,n2,n3);    
            fprintf(fp, " %14.6f ", data[index] );
            }
        fprintf(fp, "\n");
        }
    else if ( dir == 3)
        {
        i1=n1/2-1;
        i2=n2/2-1;
        for (i3 = 0; i3 < n3; ++i3)
            {
            index = index_f(i1,i2,i3,n1,n2,n3);    
            fprintf(fp, " %14.6f ", data[index] );
            }
        fprintf(fp, "\n");
        }
    else
        fprintf(stderr, " Wrong value for argument 5 in plot_data_2d \n");

    fclose(fp);
}

void plot_data_2d( char* name, int n1, int n2, int n3, int n1_local,
 int n1_local_offset, int dir, double* data)
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

    /*Cerco di capire il mio propietario della
    riga in mezzo lungo la direzione 1. In questo
    modo stampo solo quella*/
    owner=npes+1;
    if ( (n1/2 > n1_local_offset) &&  (n1/2 <= n1_local_offset + n1_local) ) {
        owner = mype;
    }

    fp = fopen (buf, "w");
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
    else if ( dir == 2)
        {
        i2=n2/2-1;
        /*Creo i buffer per le varie slice*/
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
            fp = fopen (buf, "w"); /*Scrivo i dati effetivamente sul file*/
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
        }
    else if ( dir == 3)
        {
        i3=n3/2-1;
        /*Creo i buffer per le varie slice*/
        sizes = (int*)malloc(npes*sizeof(int));
        displ = (int*)calloc(npes,sizeof(int));

        buffer= (double*)malloc(n1*n2*sizeof(double));
        buffer1d= (double*)malloc(n1*sizeof(double));
        local_buffer = (double*)malloc(n1_local*sizeof(double));

        MPI_Gather(&n1_local, 1, MPI_INT, sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if ( mype  ==  0) { /*Solo il master ha tutto*/
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
        /*Non so cosa a cosa serva ma l'ho commentato perche dava fastidio*/    
        //printf (" %d \n ", owner);
        if (mype == 0) { /*Scrivo effetivamente sul file*/
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
        }
    else {
        fprintf(stderr, " Wrong value for argument 5 in plot_data_2d \n");
    }
}
