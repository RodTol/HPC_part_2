FFTW_INC = /opt/area/shared/programs/Intel/fftw/3.3.10/openMPI/4.1.4/gnu/12.2.1/include
FFTW_LIB = /opt/area/shared/programs/Intel/fftw/3.3.10/openMPI/4.1.4/gnu/12.2.1/lib/
CC=gcc
CFLAGS = -O3 -Wall -I$(FFTW_INC)
LIBS = -L $(FFTW_LIB) -lfftw3 -lm
EXE = diffusion.x
SRC = diffusion.c plot_data.c fft_wrapper.c derivative.c 
HEAD = utilities.h
OBJ = $(SRC:.c=.o)


all: $(EXE)

$(EXE): $(OBJ)
	$(CC) $^ $(LIBS) -o $@

%.o : %.c $(HEAD)
	$(CC) $(CFLAGS) -c $<

flush:
	rm -f *.dat *.png *.gif

clean: 
	rm -f *.o *.x *~
