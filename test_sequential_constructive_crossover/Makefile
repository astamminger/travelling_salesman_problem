FC = ifort

%.o: %.f90
	$(FC) -O3 -parallel -c -o $@ $<

all: scr

scr: datatypes.o scr.o
	$(FC) -O3 -parallel -o $@ datatypes.o scr.o

clean:
	rm -f *.o
	rm -f scr
	rm -f *.mod

.PHONY: all clean

