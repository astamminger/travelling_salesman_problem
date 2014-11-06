FC = ifort

%.o: %.f90
	$(FC) -O3 -parallel -c -o $@ $<

all: optimal

optimal: datatypes.o optimal.o
	$(FC) -O3 -parallel -o $@ datatypes.o optimal.o

clean:
	rm -f *.o
	rm -f optimal
	rm -f *.mod

.PHONY: all clean

