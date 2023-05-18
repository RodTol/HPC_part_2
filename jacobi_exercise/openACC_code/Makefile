CORES = 6
dim = 60
time = 10000
r = 3
c = 3

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
CFLAGS = -acc -Wall -O3 $(INC_FLAGS)
	
ifdef flags	
	ifeq ($(flags), debug)
		CFLAGS+= -DDEBUG
	endif
endif

all: $(EXE)

$(EXE): $(OBJS)
	$(CC) $(LINKFLAGS) -o $(EXE) $(OBJS)

$(BUILD_DIR)/%.c.o: %.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CFLAGS) -c $< -o $@

run: clean all
	mpirun -np $(CORES) $(EXE) $(dim) $(time) $(r) $(c)

.PHONY: clean

clean:
	rm -f $(OBJS) $(EXE)
	
flush:
	rm -f *.dat

plot:
	@gnuplot -p plot.plt

MKDIR_P ?= mkdir -p
