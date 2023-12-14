#if using f2py
#f2py -c neigh_finder.f90 tct.f90 topo_analysis.f90 -m satk
.SUFFIXES:
OBJS = nf.o topo_analysis.o precision.o io.o tct.o
FORTRAN = gfortran
CFLAGS = -c

polikit: $(OBJS) polikit.f90
	$(FORTRAN) $^ $(OBJS) -o polikit -fcheck='all'
precision.o:precision.f90
	$(FORTRAN) $^ $(CFLAGS)
io.o:io.f90
	$(FORTRAN) $^ $(CFLAGS)
nf.o:nf.f90
	$(FORTRAN) $^ $(CFLAGS)
%.o:%.f90
	$(FORTRAN) $^ $(CFLAGS)
# %.o:%.f90
# 	$(FORTRAN) $^ $(CFLAGS)
.PHONY: clean
clean:
	rm *.o *.mod
