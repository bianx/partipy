.POSIX:
.SUFFIXES: .o .c
.PHONY: clean lint

include ../conf.mk
M = \
main\

M_CFLAGS = -I$(PREFIX)/include
M_LDFLAGS = -L$(PREFIX)/lib -lbarnes-hut -lm

all: $M
main.o: double/real.h
main.o: float/real.h

.c.o:; $(CC) -c $(CFLAGS) $(M_CFLAGS) $<
.o:; $(LINK) $< $(LDFLAGS) $(M_LDFLAGS) -o $@
.c:
%: %.c

lint:; make CFLAGS='-Wall -Wextra -g -O2'

clean:; rm -f -- $M $(M:=.o)
