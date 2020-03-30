.POSIX:
.SUFFIXES: .o .c
.PHONY: clean lint

CC = c99
CFLAGS = -g -Ofast
LINK = $(CC)
M = \
main\
main2\

GSL_CFLAGS = `pkg-config --cflags gsl`
GSL_LDFLAGS = `pkg-config --libs gsl`
M_LDFLAGS = -lm

all: $M
main2.o: double/real.h
main2.o: float/real.h

.c.o:; $(CC) -c $(CFLAGS) $(GSL_CFLAGS) $<
.o:; $(LINK) $< $(LDFLAGS) $(GSL_LDFLAGS) $(M_LDFLAGS) -o $@
.c:
%: %.c

lint:; make CFLAGS='-Wall -Wextra -g -O2'

clean:; rm -f -- $M
