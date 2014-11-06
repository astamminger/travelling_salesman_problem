FC = ifort
#FC = gfortran

%.o: %.f90
	$(FC) -O3 -openmp -c -o $@ $<

all: travelling_salesman gencity

travelling_salesman: datatypes.o travelling_salesman.o
	$(FC) -O3 -openmp -o $@ datatypes.o travelling_salesman.o

gencity: gencity.o
	$(FC) -O3 -openmp -o gencity gencity.o

clean:
	rm -f *.o
	rm -f travelling_salesman gencity
	rm -f *.mod

.PHONY: all clean

