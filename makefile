# Image: libPolynomial.A.dylib 1.0

# Common options
CFLAGS_COM = -O3 -Wall -std=gnu99

# Link against default libm
CFLAGS_LIBM  = -fPIC $(CFLAGS_COM) $(CPPFLAGS)
LDFLAGS_LIBM = -lm

# Link against Intel MKL
MKLROOT = /opt/intel/mkl
CFLAGS_MKL  = $(CFLAGS_COM) -m64 -fPIC -I${MKLROOT}/include
LDFLAGS_MKL = ${MKLROOT}/lib/libmkl_intel_lp64.a ${MKLROOT}/lib/libmkl_sequential.a ${MKLROOT}/lib/libmkl_core.a -lpthread -lm -ldl

# LMATH = LIBM
LMATH = MKL

MOD_CFLAGS = $(CFLAGS_$(LMATH))
MOD_LDFLAGS = $(LDFLAGS_$(LMATH))

all: dylib tests

dylib: Polynomial.c Polynomial.h
	clang $(MOD_CFLAGS) -dynamiclib Polynomial.c -current_version 1.0 -compatibility_version 1.0 -fvisibility=hidden -o libPolynomial.A.dylib $(MOD_LDFLAGS)
	ln -sf ./libPolynomial.A.dylib ./libPolynomial.dylib

tests: tests.c
	clang $(MOD_CFLAGS) tests.c libPolynomial.A.dylib -o tests $(MOD_LDFLAGS)

install:
	cp -f Polynomial.h /usr/local/include/
	cp -f libPolynomial.A.dylib /usr/local/lib/
	ln -sf /usr/local/lib/libPolynomial.A.dylib /usr/local/lib/libPolynomial.dylib

uninstall:
	rm -f /usr/local/lib/libPolynomial.A.dylib
	rm -f /usr/local/lib/libPolynomial.dylib
	rm -f /usr/local/include/Polynomial.h

clean:
	rm -f libPolynomial.dylib libPolynomial.A.dylib tests
