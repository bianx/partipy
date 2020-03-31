.POSIX:
.SUFFIXES: .o .c
.PHONY: clean lint

include ../conf.mk
M = \
main\

M_LDFLAGS = -lm

all: $M
main.o: double/real.h
main.o: float/real.h

.c.o:; $(CC) -c $(CFLAGS) $<
.o:; $(LINK) $< $(LDFLAGS) $(M_LDFLAGS) -o $@
.c:
%: %.c

lint:; make CFLAGS='-Wall -Wextra -g -O2'

clean:; rm -f -- $M $(M:=.o)
