.POSIX:
.PHONY: clean run all

CC99 = mpicc

M = vortex
Q_FLAGS = -lm -O3

all: $M
.c:
	. ba.setenv && \
	CC99="$(CC99)" qcc $(CFLAGS) $(LDFLAGS) $(Q_FLAGS) $< -o $@

run: all
	./vortex
clean:; rm -f $M


