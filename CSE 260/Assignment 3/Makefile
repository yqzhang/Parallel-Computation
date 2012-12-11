HOST = $(shell hostname)
BANG   =  $(shell expr match `hostname` ccom-bang)
BANG-COMPUTE   =  $(shell expr match `hostname` compute)
LILLIPUT   =  $(shell expr match `hostname` lilliput)



ifneq ($(BANG), 0)
PUB     = /share/class/public/cse260-fa12
include $(PUB)/Arch/arch.gnu.generic
else
ifneq ($(BANG-COMPUTE), 0)
PUB     = /share/class/public/cse260-fa12
include $(PUB)/Arch/arch.gnu.generic
else
ifneq ($(LILLIPUT), 0)
PUB	= /class/public/cse260-fa12
include $(PUB)/Arch/arch.intel.generic
else
# PUB = /Users/baden/lib
include $(PUB)/Arch/arch.gnu
# include $(PUB)/Arch/arch.gnu-4.5
endif
endif
endif
#
# Add symbol table information for gdb/cachegrind
ifeq ($(debug), 1)
        CFLAGS += -g
        LDFLAGS += -g
        C++FLAGS += -g
endif   


# Add symbol table information for gprof
ifeq ($(gprof), 1)
        CFLAGS += -g -pg
        C++FLAGS += -g -pg
        LDFLAGS += -g -pg
endif

# If you want to compile for single precision,
# specify single=1 on the "make" command line
ifeq ($(single), 1)
else
    C++FLAGS += -D_DOUBLE
    CFLAGS += -D_DOUBLE
endif


# If you want to compile so that you call the plotter for
# flattened 2D arrays (implemented as 1D arrays)
# specify flattened=1 on the "make" command line
ifeq ($(flattened), 1)
    C++FLAGS += -DPLOT1D
    CFLAGS += -DPLOT1D
endif

# If you want to use restrict pointers, make restrict=1
# This applies to the hand code version
ifeq ($(restrict), 1)
    C++FLAGS += -D__RESTRICT
    CFLAGS += -D__RESTRICT
ifneq ($(CARVER), 0)
    C++FLAGS += -restrict
    CFLAGS += -restrict
endif
endif


#DEBUG += -DDEBUG


app:		apf

OBJECTS = apf.o solve.o splot.o cmdLine.o Report.o printMat.o
ifeq ($(no-mpi),1)
OBJECTS += Timer.o
endif

apf:	        $(OBJECTS) 
		$(C++LINK) $(LDFLAGS) -o $@ $(OBJECTS)  $(LDLIBS)

clean:	
	$(RM) *.o apf;
	$(RM) core;
