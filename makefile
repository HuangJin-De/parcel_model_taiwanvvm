# makefile 
FC = mpiifort
CHECK = -CB -g -traceback -check all,noarg_temp_created -debug all
FCFLAGS = -r8 -O3 -free -mcmodel=large -heap-arrays 10 -shared-intel -fp-model precise -qmkl
ifdef DEBUG
  FCFLAGS += $(CHECK)
endif
FINCLUDE = -I/opt/libs-intel-oneapi/netcdf-4.7.4/include
LDLIBS = -L/opt/libs-intel-oneapi/netcdf-4.7.4/lib -lnetcdff 


## code paths
#VPATH = ./src/
#
## objects
#
#LIST = cal_pblh.F cal_cape.F
#a.out: cal_pblh.o cal_cape.o
#cal_pblh.o: cal_pblh.F cal_cape.o
#cal_cape.o: cal_cape.F 
#
#LIST_o = $(LIST:.f=.o)
#target = a.out 

# code paths
VPATH = ./

# objects

#LIST = cal_pblh.f cal_cape.f
#a.out: cal_pblh.o cal_cape.o
#cal_pblh.o: cal_pblh.f cal_cape.o
#cal_cape.o: cal_cape.f


#LIST = cal_mean.f
#a.out: cal_mean.o
#cal_mean.o: cal_mean.f

LIST = cal_wpdf.f
a.out: cal_wpdf.o
cal_mean.o: cal_wpdf.f

LIST_o = $(LIST:.f=.o)
target = a.out 

all: $(target)

#$(LIST_o): %.o: %.F
#	$(FC) $(FCFLAGS) $(FINCLUDE) -c $<

$(LIST_o): %.o: %.f
	$(FC) $(FCFLAGS) $(FINCLUDE) -c $<

$(target) : $(LIST_o)
	$(FC) $(FCFLAGS) $(FINCLUDE) $^ -o $@ $(LDLIBS)

clean:
	rm -rf *.o *.mod a.out


