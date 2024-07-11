# still to-do
# compute bgcorr vs. sky selection

ALL = ../bin/wfpc2distort ../bin/wfpc2fitdistort ../bin/wfpc2mask wfpc2phot.o ../bin/wfpc2showpsf wfpc2filters.o wfpc2distort.o ../bin/wfpc2unmask
#ALL = ../bin/wfpc2fitdistort ../bin/wfpc2mask wfpc2phot.o ../bin/wfpc2showpsf ../bin/wfpc2fakelist wfpc2filters.o wfpc2distort.o ../bin/wfpc2undistort

ifdef USETT
ALL += ../bin/wfpc2ttpar ../bin/wfpc2makepsf
endif

all: $(ALL)

wfpc2phot.o: wfpc2phot.c wfpc2psfdata.h ../include/fits.h wfpc2distort.h wfpc2filters.h ../dolphot_defs.h
	$(CC) $(CFLAGS) -c -o $@ $<

../bin/wfpc2distort: wfpc2distort_main.c wfpc2distort.o wfpc2distort.h wfpc2filters.o wfpc2filters.h ../include/dolphot.h
	$(CC) $(CFLAGS) -o $@ $< ../dolphot_lib.o wfpc2distort.o wfpc2filters.o $(LIBS)

../bin/wfpc2fitdistort: wfpc2fitdistort.c wfpc2distort.o wfpc2distort.h wfpc2filters.o wfpc2filters.h ../include/dolphot.h
	$(CC) $(CFLAGS) -o $@ $< ../dolphot_lib.o wfpc2distort.o wfpc2filters.o $(LIBS)

../bin/wfpc2undistort: wfpc2undistort.c wfpc2distort.o wfpc2distort.h wfpc2filters.o wfpc2filters.h ../include/dolphot.h
	$(CC) $(CFLAGS) -o $@ $< ../dolphot_lib.o wfpc2distort.o wfpc2filters.o $(LIBS)

../bin/wfpc2mask: wfpc2mask.c ../dolphot_lib.o ../fits_lib.o ../include/fits.h ../include/dolphot.h wfpc2psfdata.h wfpc2distort.o
	$(CC) $(CFLAGS) -o $@ $< ../dolphot_lib.o ../fits_lib.o wfpc2distort.o $(LIBS)

ifdef USETT
../bin/wfpc2ttpar: wfpc2ttpar.c
	$(CC) $(CFLAGS) -o $@ $<

../bin/wfpc2makepsf: wfpc2makepsf.c ../dolphot_lib.o ../fits_lib.o wfpc2psfdata.h ../include/fits.h ../include/dolphot.h wfpc2distort.o
	$(CC) $(CFLAGS) -o $@ $< ../dolphot_lib.o ../fits_lib.o wfpc2distort.o $(LIBS)
endif

../bin/wfpc2showpsf: wfpc2showpsf.c ../dolphot_lib.o wfpc2psfdata.h ../include/dolphot.h
	$(CC) $(CFLAGS) -o $@ $< ../dolphot_lib.o $(LIBS)

wfpc2filters.o: wfpc2filters.c
	$(CC) $(CFLAGS) -c -o $@ $<

wfpc2distort.o: wfpc2distort.c
	$(CC) $(CFLAGS) -c -o $@ $<

../bin/wfpc2fakelist: wfpc2fakelist.c wfpc2filters.h wfpc2filters.o
	$(CC) $(CFLAGS) -o $@ $< wfpc2filters.o $(LIBS)

../bin/wfpc2unmask: wfpc2unmask.c ../dolphot_lib.o ../fits_lib.o ../include/fits.h ../include/dolphot.h wfpc2psfdata.h wfpc2distort.o
	$(CC) $(CFLAGS) -o $@ $< ../dolphot_lib.o ../fits_lib.o wfpc2distort.o $(LIBS)

clean:
	rm $(ALL)

