# Makefile 
#

CFITSIOINC := -I$(ASTROSOFT)/include
CFITSIOLIB := -L$(ASTROSOFT)/lib
CC = gcc 
CFLAGS = -fPIC 
LIBS = -lcfitsio 

all:
	$(CC) ortfits.c -o ortfits $(CFITSIOINC) $(CFITSIOLIB) $(CFLAGS) $(LIBS)
	$(CC) gmrtfits.c -o gmrtfits $(CFITSIOINC) $(CFITSIOLIB) $(CFLAGS) $(LIBS)

ortfits:
	$(CC) ortfits.c -o ortfits $(CFITSIOINC) $(CFITSIOLIB) $(CFLAGS) $(LIBS)
	
gmrtfits:
	$(CC) gmrtfits.c -o gmrtfits $(CFITSIOINC) $(CFITSIOLIB) $(CFLAGS) $(LIBS)
clean:
	-rm -f ortfits gmrtfits 

