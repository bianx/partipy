.POSIX:
.SUFFIXES: .o .c
.PHONY: clean lint install
include ../../conf.mk

O = \
main.o\

L = libbarnes-hut.a
H = barnes-hut.h

install: $L $H
	mkdir -p -- "$(PREFIX)/lib" "$(PREFIX)/include"
	cp -- $H "$(PREFIX)/include"
	cp -- $L "$(PREFIX)/lib"
$L: $O; ar rv $@ $O && ranlib $@

.c.o:; $(CC) -c $(CFLAGS) $< -o $@
main.o: barnes-hut.h
include dep.mk

lint:; make CFLAGS='-Wall -Wextra -g -O2'
dep:; for i in *.c; do $(CPP) -MM $$i; done > dep.mk

clean:; rm -f -- $O $L
