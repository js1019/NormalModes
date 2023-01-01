
include ../makefile.in

exe  = ../bin/plmvcg_$(PLAT).out

sourcef90 = parmetis_interface.f90 lapack_interface.f90 mod_string.f90 mod_para.f90\
mod_datatype.f90 mod_utility.f90 mod_cg_datatype.f90 mod_geometry.f90\
mod_cg_models.f90 mod_cg_create_matrix.f90\
mod_exafmm.f90 mod_matvec.f90 mod_pevsl.f90 mainnm.f90

objectf90 = $(sourcef90:.f90=.o)

c_obj = parmetis_header.o

default:$(c_obj) $(objectf90) $(iterobj) $(exe) 

$(exe): $(sourcef90) $(objectf90) 
	$(PFC) $(OPTL) $(PFFLAGS) -o $@ $(c_obj) $(objectf90) $(LORDERINGS) $(LIBEVSL) $(iterobj) $(LDFLAGS) $(LIBFMM)


%.o: %.f90 
	$(PFC) $(FDEBUG) $(OPTL) -o $@ $<

parmetis_header.o: parmetis_header.c
	$(PCC) -c $(IPARMETIS) $(IMETIS) $<



clean: 
	$(RM) -f *.o *.mod ../bin/plmvcg_$(PLAT).out
