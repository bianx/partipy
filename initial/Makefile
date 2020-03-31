.POSIX:
.SUFFIXES: .o .c
.PHONY: clean lint
include ../conf.mk

M = \
grid\

M_CFLAGS = -I../lib/vof
M_LDFLAGS = -L../lib/vofi -L../lib/vof -lvof -lvofi -lm

all: $M
.c.o:; $(CC) -c $(M_CFLAGS) $(CFLAGS) $<
.o:; $(LINK) $< $(M_LDFLAGS) $(LDFLAGS) -o $@
.c:
%: %.c

lint:; make CFLAGS='-Wall -Wextra -g -O2'

clean:; rm -f -- $M $(M:=.o)
