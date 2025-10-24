#*****************************************************
#
#                       SOLVER OPTIONS
#
#*****************************************************

#Name of cluster: Peter, rocket
CLUSTER =       Peter


#*****************************************************
#
#                       COMPILER OPTIONS
#
#*****************************************************
NAME    =       NS_solver

# Local linux machine
ifeq ($(CLUSTER),Peter)
COMP    =       gcc
FLAGDEBUG       =       -pedantic -Wall -fcheck=all -fbacktrace -g3 -fsanitize=address
FLAGS   =       -O3 -ffast-math
FLAG_OMP=       -fopenmp
IDIR    =       -I./include
LFFTW   =       -lfftw3_threads -lfftw3
FFTWLDIR=       -L/usr/lib -L/usr/local/lib
LIBS    =       -lm -lpthread
endif
#

# Rocket cluster
ifeq ($(CLUSTER),rocket)
COMP    =       gcc
FLAGDEBUG       =       -pedantic -Wall -fcheck=all -fbacktrace -g3 -fsanitize=address
FLAGS   =       -O3 -ftree-vectorize -ffast-math
FLAG_OMP=       -fopenmp
IDIR    =       -I./include -I/usr/local/include/ -I/usr/include/ -I/mnt/storage/apps/eb/software/FFTW/3.3.10-gompi-2021b/include
LFFTW   =       -lfftw3_threads -lfftw3
FFTWLDIR=       -L/usr/lib/ -L/usr/local/lib/ -L/mnt/storage/apps/eb/software/FFTW/3.3.10-gompi-2021b/lib
LIBS    =       -lm -lpthread
endif

# -- Source & object files
SRCS = $(NAME).c params.c wavenumbers.c IC.c fields.c field_ops.c  NS_drive.c outputs.c
OBJS = $(SRCS:.c=.o)

# Compiler flags for all builds
CFLAGS = $(FLAGS) $(FLAG_OMP) $(IDIR)

# Linker flags
LDFLAGS = $(FLAGS) $(FLAG_OMP) $(FFTWLDIR) $(LFFTW) $(LIBS)

# -- Default target ---
all: $(NAME)

# --- Linking step --- CORRECTED ORDER
$(NAME): $(OBJS)
	$(COMP) $(OBJS) -o $(NAME) $(LDFLAGS)

# --- Compilation rule for .c → .o --- CORRECTED
%.o: %.c
	$(COMP) $(CFLAGS) -c $< -o $@

# --- Debug target ---
debug: CFLAGS = $(FLAGDEBUG) $(FLAG_OMP) $(IDIR)
debug: LDFLAGS = $(FLAGDEBUG) $(FLAG_OMP) $(FFTWLDIR) $(LFFTW) $(LIBS)
debug: $(NAME)_debug

$(NAME)_debug: $(OBJS)
	$(COMP) $(OBJS) -o $(NAME)_debug $(LDFLAGS)

# --- Clean up ---
clean:
	rm -f *.o $(NAME) $(NAME)_debug

# --- Test FFTW threading ---
test_fftw: test_fftw_threads.c
	$(COMP) test_fftw_threads.c -o test_fftw $(CFLAGS) $(LDFLAGS)
	./test_fftw
