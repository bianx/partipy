.POSIX:
.SUFFIXES: .o .c
.PHONY: clean lint

CC = c99
CFLAGS = -g -Og
LINK = $(CC)
M = \
main\

GSL_CFLAGS = `pkg-config --cflags gsl`
GSL_LDFLAGS = `pkg-config --libs gsl`
M_LDFLAGS = -lm

all: $M
.c.o:; $(CC) -c $(CFLAGS) $(GSL_CFLAGS) $<
.o:; $(LINK) $< $(LDFLAGS) $(GSL_LDFLAGS) $(M_LDFLAGS) -o $@
.c:
%: %.c

lint:; make CFLAGS='-Wall -Wextra -g -O2'

clean:; rm -f -- $M
