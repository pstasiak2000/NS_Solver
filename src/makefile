# --- Configuration --- 
NAME	=	NS_solver
CFLAGS	=	-O2 -Wall -I../include
LDFLAGS	=	-lfftw3 -lm

# -- Source & object files
SRCS = $(NAME).c wavenumbers.c IC.c fields.c vector_ops.c
OBJS = $(SRCS:.c=.o)

# -- Default target ---
all: $(NAME)

# --- Linking step ---
$(NAME): $(OBJS)
	gcc $(CFLAGS) $(OBJS) -o $(NAME) $(LDFLAGS)

# --- Compilation rule for .c → .o ---
%.o: %.c
	gcc $(CFLAGS) -c $< -o $@

# --- Debug target ---
debug: CFLAGS += -g -fsanitize=address -fno-omit-frame-pointer
debug: $(NAME)_debug

$(NAME)_debug: $(OBJS)
	gcc $(CFLAGS) $(OBJS) -o $(NAME)_debug $(LDFLAGS)

# --- Clean up ---
clean:
	rm -f $(OBJS) $(NAME) $(NAME)_debugF