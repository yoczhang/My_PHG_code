FILES =	Makefile test5.c

default: lib test5

all: lib test5

examples.zip: ${FILES}
	@zip -9 -u -y $@ $^

clean:
	-/bin/rm -f *.o *.tmp plot.eps PI[0-9]* \
		*.G.* *.A.* *.Aalpha.* *.Abeta.* \
		*.Gx.* *.Gy.* *.Gz.* *.x.* *.y.* *.z.* *.b.* *.x0.* \
		test5

distclean: clean
	-/bin/rm -f test1 test3 test4 test5 matrix_test mesh*.* 

lib:
	@(cd ../../src; $(MAKE))

include ../../Makefile.inc

#test1.o:		../../src/libphg${LIB_SUFFIX} test1.c
#matrix_test.o:		../../src/libphg${LIB_SUFFIX} matrix_test.c
#test3.o:		../../src/libphg${LIB_SUFFIX} test3.c
#legendre_rule.o:        ../../src/libphg${LIB_SUFFIX} legendre_rule.c
#test4.o:		../../src/libphg${LIB_SUFFIX} test4.c
test5.o:		../../src/libphg${LIB_SUFFIX} test5.c
.PHONY: default all clean distclean lib
