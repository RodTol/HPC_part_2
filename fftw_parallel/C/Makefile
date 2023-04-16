FFTW_INC = /usr/local/include
FFTW_LIB = /usr/local/lib
CC=mpicc
CFLAGS = -O3 -Wall -I$(FFTW_INC)
LIBS = -L$(FFTW_LIB) -lfftw3_mpi -lfftw3 -lm
EXE = diffusion.x
SRC = diffusion.c plot_data.c fft_wrapper.c derivative.c 
HEAD = utilities.h
OBJ = $(SRC:.c=.o)
CORES = 4

all: $(EXE)

$(EXE): $(OBJ)
	$(CC) $^ $(LIBS) -o $@

%.o : %.c $(HEAD)
	$(CC) $(CFLAGS) -c $<

run:
	mpirun -np $(CORES) diffusion.x

plot:
	gnuplot animate.plt

flush:
	rm -f *.dat
	rm -f *.gif
	rm -f *.png
reset:
	rm -f *.out

clean: 
	rm -f *.o *.x *~
