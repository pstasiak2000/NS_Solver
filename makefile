#*****************************************************
#
#			SOLVER OPTIONS
#
#*****************************************************

#Name of cluster: Peter, rocket
CLUSTER	=	Peter


#*****************************************************
#
#			COMPILER OPTIONS
#
#*****************************************************
NAME	=	NS_solver

# Local linux machine
ifeq ($(CLUSTER),Peter)
COMP	=	gcc
FLAGDEBUG	=	-pedantic -Wall -fcheck=all -fbacktrace -g3  
FLAGS	=	-O3 -ftree-vectorize -ffast-math -lpthread -lm
FLAG_OMP=	-fopenmp -lfftw3_threads
IDIR	=	-I./include
LFFTW	=	-lfftw3 -lfftw3_threads
FFTWLDIR=	-L/usr/lib -L/usr/local/lib
endif
#

# Rocket cluster
ifeq ($(CLUSTER),rocket)
COMP	=	gcc
FLAGDEBUG	=	-pedantic -Wall -fcheck=all -fbacktrace -g3
FLAGS	=	-O3 -ftree-vectorize -ffast-math -lpthread -lm
FLAG_OMP=	-fopenmp -lfftw3_threads
IDIR	=	-I./include -I/usr/local/include/ -I/usr/include/ -I/mnt/storage/apps/eb/software/FFTW/3.3.10-gompi-2021b/include
LFFTW	=	-lfftw3 -lfftw3_threads
FFTWLDIR=	-L/usr/lib/ -L/usr/local/lib/ -L/mnt/storage/apps/eb/software/FFTW/3.3.10-gompi-2021b/lib
endif

# -- Source & object files
SRCS = $(NAME).c params.c wavenumbers.c IC.c fields.c field_ops.c  NS_drive.c outputs.c
OBJS = $(SRCS:.c=.o)

# -- Default target ---
all: $(NAME)

# --- Linking step ---
$(NAME): $(OBJS)
	$(COMP) $(CFLAGS) $(OBJS) -o $(NAME) $(FLAGS) $(FLAG_OMP) $(FFTWLDIR) $(IDIR) $(LFFTW)

# --- Compilation rule for .c → .o ---
%.o: %.c
	$(COMP) $(CFLAGS) -c $< -o $@ $(FLAGS) $(FLAG_OMP) $(FFTWLDIR) $(IDIR) $(LFFTW)

# --- Debug target ---
debug: CFLAGS += -g -fsanitize=address -fno-omit-frame-pointer
debug: $(NAME)_debug

$(NAME)_debug: $(OBJS)
	gcc $(FLAGDEBUG) $(OBJS) -o $(NAME)_debug $(FLAGS) $(FLAG_OMP) $(FFTWLDIR) $(IDIR) $(LFFTW)

# --- Clean up ---
clean:
	rm -f *.o
