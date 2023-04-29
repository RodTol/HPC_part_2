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
CFLAGS = -Wall -O3 $(INC_FLAGS)

#extra flags
LINKFLAGS = 
OPENBLAS = -I ${OPENBLAS_HOME}/include/ -L ${OPENBLAS_HOME}/lib -lopenblas -lgfortran
CUBLAS = -lcublas

#make [command] flags=dgemm for example
ifdef flags
	ifeq ($(flags), dgemm)
		CFLAGS += -DDGEMM
		LINKFLAGS += $(OPENBLAS)
	else ifeq ($(flags), gpu)
		SRCS += $(SRC_DIRS)/gpu_computation.cu
		OBJS += $(BUILD_DIR)/gpu_computation.o
		CFLAGS += -DGPU
		LINKFLAGS += $(CUBLAS)
	endif
	
	ifeq ($(flags), debug)
		CFLAGS+= -DDEBUG
	endif
endif

all: $(EXE)

$(EXE): $(OBJS)
	$(CC) -o $(EXE) $(OBJS)

$(BUILD_DIR)/%.c.o: %.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CFLAGS) $(LINKFLAGS) -c $< -o $@

run: clean all
	mpirun -np $(CORES) $(EXE) $(N)

.PHONY: clean

clean:
	rm -f $(OBJS) $(EXE)
	rm -r $(BUILD_DIR)

flush:
	rm -f *.dat

MKDIR_P ?= mkdir -p

