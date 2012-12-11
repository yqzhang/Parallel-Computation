# Makefile
HOST = $(shell hostname)
BANG = $(shell expr match `hostname` ccom-bang)

ifneq ($(BANG), 0)
PUB	= /share/class/public/cse260-fa12
include $(PUB)/Arch/arch.gnu.generic
else
endif

ifeq ($(single), 1)
else
	C++FLAGS += -D_DOUBLE
	CFLAGS += -D_DOUBLE
endif

ifeq ($(debug), 1)
	C++FLAGS += -ggdb
	CFLAGS += -ggdb
	LDFLAGS += -ggdb
endif

app:	driver

OBJECTS = driver.o fft_mpi.o fft_serial.o

driver:	$(OBJECTS)
	$(CLINK) $(LDFLAGS) -o $@ $(OBJECTS) $(LDLIBS)

clean:
	$(RM) *.o driver;
