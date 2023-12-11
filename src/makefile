#if using f2py
#f2py -c neigh_finder.f90 tct.f90 topo_analysis.f90 -m satk
.SUFFIXES:
OBJS = neigh_finder.o topo_analysis.o precision.o io.o tct.o
FORTRAN = gfortran
CFLAGS = -c

a.exe: $(OBJS) main.f90
	$(FORTRAN) main.f90 $(OBJS) -o a.exe -fcheck='all'
%.o:%.f90
	$(FORTRAN) $^ $(CFLAGS)
.PHONY: clean
clean:
	rm *.o *.mod
