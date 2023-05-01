CORES = 6
N = 6

#target
BUILD_DIR = ./build
EXE = a.out

SRC_DIRS ?= ./src

#sources
SRCS = $(shell find $(SRC_DIRS) -name *.c)
OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)

#headers
INC_DIRS = $(shell find $(SRC_DIRS) -type d)
INC_FLAGS = $(addprefix -I,$(INC_DIRS))

#compiler and standard flags
CC = mpicc
NVCC = nvcc
CFLAGS = -Wall -O3 $(INC_FLAGS)

#extra flags
LINKFLAGS = 
LOPENBLAS = -I ${OPENBLAS_HOME}/include/ -L ${OPENBLAS_HOME}/lib -lopenblas -lcblas

#for compilation on local system
CUDA_HOME = /usr/lib/cuda
LCUBLAS = -L${CUDA_HOME}/lib64/ -lcudart -lcublas

#make [command] compilation=dgemm for example
ifdef compilation
	ifeq ($(compilation), dgemm)
		CFLAGS += -DDGEMM
		LINKFLAGS += $(LOPENBLAS)
	else ifeq ($(compilation), gpu)
		OBJS += $(BUILD_DIR)/gpu_computation.c.o
		CFLAGS += -DGPU
		LDLIBS += $(LCUBLAS)
	endif
endif
	
ifdef flags	
	ifeq ($(flags), debug)
		CFLAGS+= -DDEBUG
	endif
endif

all: $(EXE)

$(EXE): $(OBJS)
	$(CC) $(CFLAGS) $(LINKFLAGS) -o $(EXE) $(OBJS) $(LDLIBS) 
 
$(BUILD_DIR)/%.c.o: %.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CFLAGS) $(LINKFLAGS) -c $< -o $@

$(BUILD_DIR)/gpu_computation.c.o : $(SRC_DIRS)/gpu_computation.cu
	$(NVCC) $ -c $< -o $@


run: clean all
	mpirun -np $(CORES) $(EXE) $(N)

.PHONY: clean

clean:
	rm -f $(OBJS) $(EXE)
	
flush:
	rm -f *.dat

MKDIR_P ?= mkdir -p

