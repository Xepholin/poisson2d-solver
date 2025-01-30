CC=gcc
CFLAGS=-g -Wall -Wextra -march=native

OFLAGS=-O3

INTEL_LFLAGS=-qmkl
GNU_LFLAGS=-lm -fopenmp

FILES=main.c

all: poisson

run: poisson
	./poisson 4 150 0.0000001

poisson: $(FILES)
ifeq ($(CC),icc)
	$(CC) $(CFLAGS) $(OFLAGS) $(FILES) -o $@ $(INTEL_LFLAGS)
else
ifeq ($(CC),icx)
	$(CC) $(CFLAGS) $(OFLAGS) $(FILES) -o $@ $(INTEL_LFLAGS)
else
ifeq ($(CC),gcc)
	$(CC) $(CFLAGS) $(OFLAGS) $(FILES) -o $@ $(GNU_LFLAGS)
else
ifeq ($(CC),clang)
	$(CC) $(CFLAGS) $(OFLAGS) $(FILES) -o $@ $(GNU_LFLAGS)
else
ifeq ($(CC),aocc)
	$(CC) $(CFLAGS) $(OFLAGS) $(FILES) -o $@ $(GNU_LFLAGS)
else
	@echo "ERROR: no compiler specified using CC. Possible values for CC: gcc, icc, icx, clang"
endif
endif
endif
endif
endif

clean:
	@rm -Rf poisson

