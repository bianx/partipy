.POSIX:
.SUFFIXES: .o .c
.PHONY: clean lint
include ../../conf.mk

O = \
main.o\

L = libtree.a
H = tree.h

$L: $O; ar rv $@ $O && ranlib $@

.c.o:; $(CC) -c $(CFLAGS) $< -o $@
include dep.mk

lint:; make CFLAGS='-Wall -Wextra -g -O2'
dep:; for i in *.c; do $(CPP) -MM $$i; done > dep.mk

clean:; rm -f -- $O $L
