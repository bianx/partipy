.POSIX:
.SUFFIXES: .o .c

include ../../conf.mk
include obj.mk
L = libvofi.a
H = vofi.h

$L: $O; ar rv $@ $O && ranlib $@

.c.o:; $(CC) -c $(CFLAGS) $< -o $@
include dep.mk
lint:; make CFLAGS='-Wall -Wextra -g -O2'
dep:; for i in *.c; do $(CPP) -MM $$i; done > dep.mk
obj:; { echo O = \\; ls *.c | sed -e 's,^\.\/,,1' -e 's,\.c,.o\\,1'; echo; } > obj.mk

.PHONY: clean lint install dep obj
clean:; rm -f -- $O $L
