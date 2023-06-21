# local environment
FFTW_INC = /home/rodolfo/spack/opt/spack/linux-ubuntu22.04-skylake/gcc-12.1.0/fftw-3.3.10-tqiaj7s2nkff62gmiwmrainm5mjcedur/include
FFTW_LIB = /home/rodolfo/spack/opt/spack/linux-ubuntu22.04-skylake/gcc-12.1.0/fftw-3.3.10-tqiaj7s2nkff62gmiwmrainm5mjcedur/lib

#M100 env
#FFTW_INC = /m100_work/PROJECTS/spack/spack-0.14/install/linux-rhel7-power9le/gcc-8.4.0/fftw-3.3.8-hwlrarpm6cvjlukhfdowwveb7g7oqwgc/include
#FFTW_LIB = /m100_work/PROJECTS/spack/spack-0.14/install/linux-rhel7-power9le/gcc-8.4.0/fftw-3.3.8-hwlrarpm6cvjlukhfdowwveb7g7oqwgc/lib

n_socket=16
n_node=32

CC=mpicc
LIBS = -L$(FFTW_LIB) -lfftw3_mpi -lfftw3 -lm
EXE = a.out
SRC_DIR = src
BUILD_DIR = build
SRCS = $(shell find $(SRC_DIR) -name *.c)
OBJS = $(SRCS:%=$(BUILD_DIR)/%.o)
OBJ = $(OBJS:.c.o=.o)

INC_DIRS = $(shell find $(SRC_DIR) -type d)
INC_FLAGS = $(addprefix -I ,$(INC_DIRS))

CFLAGS = -O3 -Wall -I$(FFTW_INC) $(INC_FLAGS)
CORES = 6

ifdef flags
	ifeq ($(flags), debug)
		CFLAGS+= -DDEBUG
	endif

	ifeq ($(flags), info)
		CFLAGS+= -DPRINT_INFO
	endif
endif

all: $(EXE)

$(EXE): $(OBJ)
	$(CC) $^ $(LIBS) -o $@

$(BUILD_DIR)/%.o : %.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CFLAGS) -c $< -o $@

run: flush clean all
	mpirun -np $(CORES) -npersocket $(n_socket) -npernode $(n_node) --bind-to core $(EXE)

plot:
	gnuplot animate.plt

flush:
	rm -f data/*.dat
	rm -f images/*.gif
	rm -f images/*.png

clean: 
	rm -f $(OBJ) $(EXE)

MKDIR_P ?= mkdir -p
