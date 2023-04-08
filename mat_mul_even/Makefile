CORES = 4

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
DDGEMM = 
OPENBLAS = -I ${OPENBLAS_HOME}/include/ -L ${OPENBLAS_HOME}/lib -lopenblas -lgfortran
CUBLAS = -lcublas

ifdef flags
	ifeq ($(flags), dgemm)
        	CFLAGS += -DDGEMM
	        BLASLINK = $(OPENBLAS)
	else ifeq ($(flags), cuda)
        	CFLAGS += -DCUDA
	        BLASLINK = $(CUBLAS)
	endif
endif

all: $(EXE)

$(EXE): $(OBJS)
	$(CC) -o $(EXE) $(OBJS)

$(BUILD_DIR)/%.c.o: %.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CFLAGS) -c $< -o $@

run: clean all
	mpirun -np $(CORES) $(EXE)

.PHONY: clean

clean:
	rm -f $(OBJS) $(EXE)
	rm -r $(BUILD_DIR)

MKDIR_P ?= mkdir -p

