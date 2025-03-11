plots: rdf_data.txt rmsd_data.txt
	python3 plot.py
rdf_data.txt rmsd_data.txt : vis.exe
	./vis.exe
vis.exe : vis.f90
	gfortran vis.f90 -o vis.exe
all: plots

