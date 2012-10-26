HOST = $(shell hostname)
TRESTLES   =  $(shell expr match `hostname` trestles)
LILLIPUT   =  $(shell expr match `hostname` lilliput)
CSECLASS   =  $(shell expr match `hostname` cseclass)
BANG   =  $(shell expr match `hostname` ccom-bang)


ifneq ($(LILLIPUT), 0)
PUB	= /class/public/cse260-fa12
include $(PUB)/Arch/arch.intel-mkl
else
ifneq ($(CSECLASS), 0)
PUB	= /class/public/cse260-wi12
include $(PUB)/Arch/arch.intel
else
ifneq ($(TRESTLES), 0)
PUB     = /home/baden/cse260-wi12
include $(PUB)/Arch/arch.intel
else
ifneq ($(BANG), 0)
PUB = /share/class/public/cse260-fa12
include $(PUB)/Arch/arch.gnu.atlas
else
PUB = /Path/To/Your/Own/machine/here
# PUB = /Users/baden/lib
include $(PUB)/Arch/arch.gnu
# include $(PUB)/Arch/arch.gnu-4.5
endif
endif
endif
endif

# XTRAFLAGS += -fsource-asm -S

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


# If you want to copy data blocks to contiguous storage
# This applies to the hand code version
ifeq ($(copy), 1)
    C++FLAGS += -DCOPY
    CFLAGS += -DCOPY
endif


# If you want to use restrict pointers, make restrict=1
# This applies to the hand code version
ifeq ($(restrict), 1)
    C++FLAGS += -D__RESTRICT
    CFLAGS += -D__RESTRICT
endif


#DEBUG += -DDEBUG


# CFLAGS += -DUSE_BLAS
# CFLAGS += -restrict

UTIL	= Timer.o
OBJ	= mmpy.o mm_unblocked.o mm_blas.o mm_blocked.o cmdLine.o verify.o $(UTIL)

mmpy:	        $(OBJ)
		$(C++LINK) $(LDFLAGS) -o $@ $(OBJ)  $(LDLIBS)



clean:	
	$(RM) *.o mmpy;
	$(RM) core;

