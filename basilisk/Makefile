.POSIX:
.PHONY: clean run all

CC99 = mpicc
CFLAGS = -Ofast -g

M = \
vortex\
wall\

Q_FLAGS = -L$(BASILISK)/gl -lglutils -lfb_glx -lGLU -lGLEW -lGL -lX11 -lm

all: $M
.c:
	CC99="$(CC99)" qcc $< $(CFLAGS) $(LDFLAGS) $(Q_FLAGS) -o $@

run: all
	./vortex
clean:; rm -f $M


