# --- Configuration --- 
NAME	=	NS_solver
COMP	=	gcc
FLAGDEBUG	=	-pedantic -Wall 
FLAGS	=	-O3 -ftree-vectorize -ffast-math -lpthread -lm
FLAG_OMP=	-fopenmp -lfftw3_threads
IDIR	=	-I./include
LFFTW	=	-lfftw3 -lfftw3_threads

# -- Source & object files
SRCS = $(NAME).c params.c wavenumbers.c IC.c fields.c field_ops.c  NS_drive.c outputs.c
OBJS = $(SRCS:.c=.o)

# -- Default target ---
all: $(NAME)

# --- Linking step ---
$(NAME): $(OBJS)
	$(COMP) $(CFLAGS) $(OBJS) -o $(NAME) $(FLAGS) $(FLAG_OMP) $(IDIR) $(LFFTW)

# --- Compilation rule for .c → .o ---
%.o: %.c
	$(COMP) $(CFLAGS) -c $< -o $@ $(FLAGS) $(FLAG_OMP) $(IDIR) $(LFFTW)

# --- Debug target ---
debug: CFLAGS += -g -fsanitize=address -fno-omit-frame-pointer
debug: $(NAME)_debug

$(NAME)_debug: $(OBJS)
	gcc $(FLAGDEBUG) $(OBJS) -o $(NAME)_debug $(FLAGS) $(FLAG_OMP) $(IDIR) $(LFFTW)

# --- Clean up ---
clean:
	rm -f *.o
