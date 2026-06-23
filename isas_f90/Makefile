NCDF_INC  = -I/opt/netcdf-ifx/include
NCDF_LIB  = -L/opt/netcdf-ifx/lib -lnetcdff -L/usr/lib/x86_64-linux-gnu -lnetcdf
LIBS      = -qmkl=sequential
LIBSFLAGS = $(NCDF_LIB) $(LIBS) -Wl,-rpath,/opt/netcdf-ifx/lib:/usr/lib/x86_64-linux-gnu


FCC   = ifx
FCOMP = $(FCC)
FLINK = $(FCC)

OPT_port1 = -i4 -r8 -extend-source
#OPT_port2 = -O0 -g -traceback -fpe0 -w

OPT_port2 = -O0 -g -traceback -w -heap-arrays 64
FOPT = $(OPT_port1) $(OPT_port2)

%.o: %.f90
	$(FCOMP) -free $(NCDF_INC) $(FOPT) -c $< -o $@

SRCTST = OA_main.f90 OA_anaarea.f90 \
	OA_ncreaddim.f90 OA_ncreadstat.f90 OA_ncreaddata.f90 \
	OA_covini.f90 OA_oversamp.f90 OA_calsol.f90 \
	OA_ncwrite_fld.f90 OA_ncwrite_dat.f90

OBJTST = $(SRCTST:.f90=.o)

OA_main: $(OBJTST)
	$(FLINK) $(OBJTST) $(LIBSFLAGS) -o $@

clean:
	rm -rf *.o $(OBJTST) OA_main
