#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>


void print_matrix(double * A, int dim_1, int dim_2, int dim_3 ) {
  int k=0;
  for(int i = 0; i < dim_1; i++ ){
    for(int j = 0; j < dim_2; j++ ){
      /*
      for (k = 0; k < dim_3; k++) {
  
      }
      */
      fprintf( stdout, "%.3g ", A[ dim_3*dim_2*i + dim_3*j + k ] );
    }
    fprintf( stdout, "\n");
  }
}

void print_linear(double* A, size_t size) {
  for (int i = 0; i < size; i++)
  {
    printf("|%.3f", A[i]);
  }
  printf("| \n");
}

void initialize(double * A, int dim_1, int dim_2, int dim_3, int offset ) {
  for(int i = 0; i < dim_1; i++ ){
    for(int j = 0; j < dim_2; j++ ){
      for (int k = 0; k < dim_3; k++) {
        A[ dim_3*dim_2*i + dim_3*j + k] =  dim_3*dim_2*i + dim_3*j + k + offset;
      }
    }
  }
}

int main (int argc, char* argv[]) {
  int n_proc_tot, irank;

  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &irank);
  MPI_Comm_size( MPI_COMM_WORLD, &n_proc_tot);

  int n1 = 4, n2=4, n3=2;
  int n1_loc = n1/n_proc_tot, n2_loc=n2/n_proc_tot;

  if( ( ( n1 % n_proc_tot ) || ( n2 % n_proc_tot ) ) && !irank ){
    fprintf( stdout, "\nN1 dimension must be multiple of the number of processes. The program will be aborted...\n\n" );
    MPI_Abort( MPI_COMM_WORLD, 1 );
  }

  size_t local_size=n1_loc*n2*n3;

  double *debug_data;
  debug_data = ( double* ) malloc( local_size * sizeof(double) );
  double *debug_columns;
  debug_columns = ( double* ) malloc( local_size * sizeof(double) );

  /*Inizializza*/
  initialize(debug_data, n1_loc,n2,n3, ((n1*n2*n3)/n_proc_tot)*irank);
  if (irank==0) printf("--- Original data ---\n");
  for (int i = 0; i < n_proc_tot; i++)
  { 
    MPI_Barrier(MPI_COMM_WORLD);
    if(irank==i) {
      sleep(0.5*irank);
      printf("I am %d : \n", irank);
      //print_matrix(debug_data, n1_loc,n2,n3);
      print_linear(debug_data, local_size);
    }
  }
  
  //Perform an Alltoall communication 
  MPI_Datatype column_block;
  MPI_Type_vector(n1_loc, n2_loc*n3, n2*n3, MPI_DOUBLE, &column_block);
  MPI_Type_commit(&column_block);

  int *send_n_of_blocks;
  int *recv_n_of_blocks;
  int *send_displacement;
  int *recv_displacement;
  MPI_Datatype * send_type;
  MPI_Datatype * recv_type;

  send_n_of_blocks = (int *)malloc(n_proc_tot * sizeof(int));
  recv_n_of_blocks = (int *)malloc(n_proc_tot * sizeof(int));
  send_displacement = (int *)malloc(n_proc_tot * sizeof(int));
  recv_displacement = (int *)malloc(n_proc_tot * sizeof(int));
  send_type = (MPI_Datatype *)malloc(n_proc_tot * sizeof(MPI_Datatype));
  recv_type = (MPI_Datatype *)malloc(n_proc_tot * sizeof(MPI_Datatype));

  for (int i = 0; i < n_proc_tot; i++)
  {
    send_n_of_blocks[i] = 1;
    recv_n_of_blocks[i] = n1_loc*n2_loc*n3;
    send_type[i] = column_block;
    recv_type[i] = MPI_DOUBLE; 
  }
  
  /*I displacement mi dicono all'i-esimo messaggio dove inizia il posto in cui 
  * metterlo o dove inizia la memoria da cui lanciare. Quindi il primo proc
  * invia a 0 e inizia a salvarsi da 0. Il secondo blocco inizia a salvarsi le cose
  * in 0+block_size ma le lancia dall'inizio del blocco, quindi n_loc*n3
  * 
  */
send_displacement[0] = 0*sizeof(double);
recv_displacement[0] = 0*sizeof(double);
for (int i = 0; i < n_proc_tot-1; i++)
{
  send_displacement[i+1] = send_displacement[i] + n2_loc*n3;//*sizeof(double);
  recv_displacement[i+1] = recv_displacement[i] + n1_loc*n2_loc*n3;//*sizeof(double); 
}

  //MPI_Alltoall(debug_data, 1, column_block, debug_columns, n1_loc*n2_loc*n3, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Alltoallv(debug_data, send_n_of_blocks, send_displacement, column_block,
   debug_columns, recv_n_of_blocks, recv_displacement, MPI_DOUBLE, MPI_COMM_WORLD);

  /*Devo usare la variante w perchè indicando i displacements con i normali indici non
  funziona, credo per il datatype che crea confusione*/
  //MPI_Alltoallw(debug_data, send_n_of_blocks, send_displacement, send_type,
   //debug_columns, recv_n_of_blocks, recv_displacement, recv_type, MPI_COMM_WORLD);
  /*Little change to see the modifications*/
  for (int i = 0; i < local_size; i++)
  {
    debug_columns[i] = debug_columns[i]+100.0;
  }
  sleep(1);
  if (irank==0) printf("--- AFTER ALL TO ALL ---\n Original matrix:\n");
  for (int i = 0; i < n_proc_tot; i++)
  { 
    MPI_Barrier(MPI_COMM_WORLD);
    if(irank==i) {
      sleep(0.5*irank);
      printf("I am %d : \n", irank);
      //print_matrix(debug_columns, n1,n2_loc,n3);
      print_linear(debug_data, local_size);
    }
  }
  
  sleep(1);
  if (irank==0) printf("--- AFTER ALL TO ALL ---\n Columns matrix:\n");
  for (int i = 0; i < n_proc_tot; i++)
  { 
    MPI_Barrier(MPI_COMM_WORLD);
    if(irank==i) {
      sleep(0.5*irank);
      printf("I am %d : \n", irank);
      //print_matrix(debug_columns, n1,n2_loc,n3);
      print_linear(debug_columns, local_size);
    }
  }




  /*Nota come facendo questo tipo di all to all, ho che i dati dopo sono sempre
  in row-major (quindi leggo i piani n2-n3) ma con n2 = n2_loc e lunghi n1*/
/*
  sleep(1);
  if (irank==0) printf("--- Matrix version of plane 0---\n");
  for (int i = 0; i < n_proc_tot; i++)
  { 
    MPI_Barrier(MPI_COMM_WORLD);
    if(irank==i) {
      sleep(0.5*irank);
      printf("I am %d : \n", irank);
      print_matrix(debug_columns, n1,n2_loc,n3);
    }
  }
*/

  MPI_Alltoallv(debug_columns, recv_n_of_blocks, recv_displacement, MPI_DOUBLE,
   debug_data, send_n_of_blocks, send_displacement, column_block, MPI_COMM_WORLD);

  /*Now i have to come back to the first order*/
  //MPI_Alltoallw(debug_columns, recv_n_of_blocks, recv_displacement, recv_type,
  // debug_data, send_n_of_blocks, send_displacement, send_type, MPI_COMM_WORLD);

  sleep(1);
  if (irank==0) printf("--- AFTER ALL TO ALL #2 ---\n Original matrix:\n");
  for (int i = 0; i < n_proc_tot; i++)
  { 
    MPI_Barrier(MPI_COMM_WORLD);
    if(irank==i) {
      sleep(0.5*irank);
      printf("I am %d : \n", irank);
      //print_matrix(debug_columns, n1,n2_loc,n3);
      print_linear(debug_data, local_size);
    }
  }

  sleep(1);
  if (irank==0) printf("--- AFTER ALL TO ALL #2 ---\n Columns matrix:\n");
  for (int i = 0; i < n_proc_tot; i++)
  { 
    MPI_Barrier(MPI_COMM_WORLD);
    if(irank==i) {
      sleep(0.5*irank);
      printf("I am %d : \n", irank);
      //print_matrix(debug_columns, n1,n2_loc,n3);
      print_linear(debug_columns, local_size);
    }
  }

  free(debug_data);
  free(debug_columns);

  MPI_Finalize();
  return 0;
}