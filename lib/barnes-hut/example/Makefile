.POSIX:
.SUFFIXES: .o .c
include ../../../conf.mk
M_CFLAGS = -I..
M_LDFLAGS = -L.. -lbarnes-hut -lm

M = \
force\
info\
print\

all: $M
.c.o:; $(CC) -c $(CFLAGS) $(M_CFLAGS) $<
.o:; $(LINK) $< $(LDFLAGS) $(M_LDFLAGS) -o $@
%: %.c

lint:; make CFLAGS='-Wall -Wextra -g -O2'

.PHONY: clean lint
clean:; rm -f -- $M
