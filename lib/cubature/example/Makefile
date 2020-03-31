.POSIX:
.SUFFIXES: .o .c
.PHONY: clean lint
include ../../../conf.mk

M = \
main\

M_CFLAGS = -I..
M_LDFLAGS = -L.. -lcubature -lm

all: $M
.c.o:; $(CC) -c $(M_CFLAGS) $(CFLAGS) $<
.o:; $(LINK) $< $(M_LDFLAGS) $(LDFLAGS) -o $@
.c:
%: %.c

lint:; make CFLAGS='-Wall -Wextra -g -O2'

clean:; rm -f -- $M $(M:=.o)
