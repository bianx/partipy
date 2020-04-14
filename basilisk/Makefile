.POSIX:
.PHONY: clean run all

CC99 = mpicc

M = vortex
M_LDFLAGS = -lm

all: $M
.c:
	. ba.setenv && \
	CC99="$(CC99)" qcc $(CFLAGS) $(LDFLAGS) $(M_LDFLAGS) $< -o $@

run: all
	./vortex
clean:; rm -f $M


