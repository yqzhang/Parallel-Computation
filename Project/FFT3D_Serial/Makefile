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

app:	driver

CFLAGS = -ggdb
OBJECTS = driver.o fft.o timer.o
CLINK = gcc
CC = gcc
LDFLAGS = -lm

driver.o: driver.c
	$(CC) $(CFLAGS) -c driver.c

fft.o: fft.c
	$(CC) $(CFLAGS) -c fft.c

timer.o: timer.c
	$(CC) $(CFLAGS) -c timer.c

driver:	$(OBJECTS)
	$(CLINK) $(LDFLAGS) -o $@ $(OBJECTS) $(LDLIBS)

clean:
	$(RM) *.o driver;
