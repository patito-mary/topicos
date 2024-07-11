CFITSIO_PREFIX = $(shell brew --prefix cfitsio)
CFLAGS += -I$(CFITSIO_PREFIX)
LDFLAGS += -L$(CFITSIO_PREFIX)/lib -lcfitsio
BASEDIR = $(shell pwd)


ALL = ../acsandersonpsf ../acsdistort ../acsfitdistort ../acsmask acsphot.o ../acsshowpsf ../acsfakelist acsfilters.o acsdistort.o ../acsundistort ../acsunmask

ifdef USETT
ALL += ../acsttpar ../acsmakepsf
endif

all: $(ALL)

acsphot.o: acsphot.c acspsfdata.h ../include/fits.h acsdistort.h acsfilters.h ../dolphot_defs.h
	$(CC) $(CFLAGS) -c -o $@ $<

../acsandersonpsf: acsandersonpsf.c ../dolphot_lib.o ../fits_lib.o acspsfdata.h ../include/fits.h ../include/dolphot.h
	$(CC) $(CFLAGS) -o $@ $< ../dolphot_lib.o ../fits_lib.o

../acsdistort: acsdistort_main.c acsdistort.o acsdistort.h acsfilters.o acsfilters.h ../include/dolphot.h
	$(CC) $(CFLAGS) -o $@ $< ../dolphot_lib.o acsdistort.o acsfilters.o $(LIBS)

../acsfitdistort: acsfitdistort.c acsdistort.o acsdistort.h acsfilters.o acsfilters.h ../include/dolphot.h
	$(CC) $(CFLAGS) -o $@ $< ../dolphot_lib.o acsdistort.o acsfilters.o $(LIBS)

../acsundistort: acsundistort.c acsdistort.o acsdistort.h acsfilters.o acsfilters.h ../include/dolphot.h
	$(CC) $(CFLAGS) -o $@ $< ../dolphot_lib.o acsdistort.o acsfilters.o $(LIBS)

../acsmask: acsmask.c ../dolphot_lib.o ../fits_lib.o acspsfdata.h ../include/fits.h ../include/dolphot.h 
	$(CC) $(CFLAGS) -o $@ $< ../dolphot_lib.o ../fits_lib.o $(LIBS)

ifdef USETT
../acsttpar: acsttpar.c
	$(CC) $(CFLAGS) -o $@ $<

../acsmakepsf: acsmakepsf.c ../dolphot_lib.o ../fits_lib.o acspsfdata.h ../include/fits.h ../include/dolphot.h
	$(CC) $(CFLAGS) -o $@ $< ../dolphot_lib.o ../fits_lib.o
endif

../acsshowpsf: acsshowpsf.c ../dolphot_lib.o acspsfdata.h ../include/dolphot.h
	$(CC) $(CFLAGS) -o $@ $< ../dolphot_lib.o $(LIBS)

acsfilters.o: acsfilters.c
	$(CC) $(CFLAGS) -c -o $@ $<

acsdistort.o: acsdistort.c
	$(CC) $(CFLAGS) -c -o $@ $<

../acsfakelist: acsfakelist.c acsfilters.h acsfilters.o
	$(CC) $(CFLAGS) -o $@ $< acsfilters.o $(LIBS)

../acsunmask: acsunmask.c ../dolphot_lib.o ../fits_lib.o ../include/fits.h ../include/dolphot.h acspsfdata.h
	$(CC) $(CFLAGS) -o $@ $< ../dolphot_lib.o ../fits_lib.o $(LIBS)

clean:
	rm $(ALL)

