#if using f2py
#f2py -c neigh_finder.f90 tct.f90 topo_analysis.f90 -m satk

a.exe: neigh_finder.o io.o topo_analysis.o precision.o
	gfortran main.f90 neigh_finder.o topo_analysis.o io.o precision.o -o a.exe
neigh_finder.o: io.o precision.o
	gfortran -c neigh_finder.f90
# tct.o: io.o precision.o:
# 	gfortran -c tct.f90
topo_analysis.o: neigh_finder.o io.o precision.o
	gfortran -c topo_analysis.f90
precision.o:
	gfortran -c precision.f90
io.o: precision.o
	gfortran -c io.f90
clean:
	rm *.o
