CORES = 4
dim = 20000
time = 1000
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
CFLAGS = -fast -acc -Minfo=accel $(INC_FLAGS)

#linkflags to libraries
LINKFLAGS = -acc -Minfo=accel

ifdef flags
	ifeq ($(flags), debug)
		CFLAGS+= -DDEBUG
	endif

	ifeq ($(flags), noprint)
		CFLAGS+= -DNOPRINT
	endif
endif

all: $(EXE)

$(EXE): $(OBJS)
	$(CC) $(LINKFLAGS) $(OBJS) -o $(EXE)

$(BUILD_DIR)/%.c.o: %.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CFLAGS) -c $< -o $@

run: clean all
	mpirun -np $(CORES) --map-by ppr:2:socket $(EXE) $(dim) $(time) $(r) $(c)
	#mpirun -np $(CORES) $(EXE) $(dim) $(time) $(r) $(c)

.PHONY: clean

clean:
	rm -f $(OBJS) $(EXE)

flush:
	rm -f initial.dat
	rm -f solution.dat
	rm -f images/*.png
	rm -f *t1
plot:
	@gnuplot -p plot.plt

MKDIR_P ?= mkdir -p

