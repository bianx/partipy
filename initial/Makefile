.POSIX:
.SUFFIXES: .o .c
.PHONY: clean lint
include ../conf.mk

M = \
grid\
ellipse\

GSL_CFLAGS = `pkg-config --cflags gsl`
GSL_LDFLAGS = `pkg-config --libs gsl`

M_CFLAGS = -I../lib/vof
M_LDFLAGS = -L../lib/vofi -L../lib/vof -lvof -lvofi -lm

all: $M
.c.o:; $(CC) -c $(M_CFLAGS) $(GSL_CFLAGS) $(CFLAGS) $<
.o:; $(LINK) $< $(M_LDFLAGS) $(GSL_LDFLAGS) $(LDFLAGS) -o $@
.c:
%: %.c

lint:; make CFLAGS='-Wall -Wextra -g -O2'

clean:; rm -f -- $M $(M:=.o)
