FILES =	Makefile test4.c

default: lib test4

all: lib test4

examples.zip: ${FILES}
	@zip -9 -u -y $@ $^

clean:
	-/bin/rm -f *.o *.tmp plot.eps PI[0-9]* \
		*.G.* *.A.* *.Aalpha.* *.Abeta.* \
		*.Gx.* *.Gy.* *.Gz.* *.x.* *.y.* *.z.* *.b.* *.x0.*

distclean: clean
	-/bin/rm -f test1 test4 matrix_test mesh*.* 
		

lib:
	@(cd ../../src; $(MAKE))

include ../../Makefile.inc

#test1.o:		../../src/libphg${LIB_SUFFIX} test1.c
#matrix_test.o:		../../src/libphg${LIB_SUFFIX} matrix_test.c
#test3.o:		../../src/libphg${LIB_SUFFIX} test3.c
#legendre_rule.o:        ../../src/libphg${LIB_SUFFIX} legendre_rule.c
test4.o:		../../src/libphg${LIB_SUFFIX} test4.c
.PHONY: default all clean distclean lib
