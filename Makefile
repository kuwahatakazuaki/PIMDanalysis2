program = run.exe
# +++ gfortran +++
fc = gfortran
fcopt =  -O2 -pipe -lblas -llapack
dirfile = /Users/kuwahatakazuaki/Program/bin/PIMDanalysis
#fcopt =  -Wall -O3 -fbacktrace -fbounds-check -lblas -llapack 
#fcopt = -Wall -O3
# +++ End gfortran +++

# +++ ifort +++
#fc = ifort
#fcopt =  -CB -traceback -fpe0
#fcopt =  -warn all -traceback
# +++ End ifort +++
objs = \
parameters.o       \
read_inp.o         \
read_coor.o        \
utility.o          \
hist1D.o           \
bond.o             \
angle.o            \
dihedral.o         \
hist2D.o           \
special_case.o      \
rotation.o         \
other_quantities.o  \
pbhpo4.o           \
periodic.o         \
dummy_atom.o       \
binary_calc.o       \
multi_bond.o       \
main.o             \
beads_expansion.o  \
#cent.o             \
#projection.o         \

module =              \
input_parameter.mod   \
utility.mod           \
calc_histogram1d.mod  \
calc_histogram2d.mod  \
mod_special_case.mod    \
mod_other_quantities.mod \
mod_periodic.mod    \
#calc_parameter.mod    \
#calc_centoroid.mod    \

%.mod : %.f90 %.o
	@true

$(program): $(objs)
	@echo
	$(fc) $(objs) -o $@ $(fcopt)
	@echo -e '\e[34m Noraml termination!!!\e[m\n'

%.o : %.f90
	@echo
	@echo ' << Compiling >>' '"'$<'"'
	$(fc) $(fcopt) -c $< -o $@


clean:
	rm -f *.o *.mod $(program)

install: $(objs)
	$(fc) $(fcopt) $(objs) -o $(program)
	cp $(program) $(dirfile)
# 	cp $(program) /Users/kuwahatakazuaki/PIMD/Analysis/Program/PIMDanalysis

