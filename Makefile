#*
#*
#* Author: Mauricio Araya Polo
#* Date: 04/2015 - present
#*
#*

INCL=.
BIN=.

CC=gnu

#For GNU compiler
ifeq ($(CC),gnu)
	CC=gcc
	CFLAGS=-g -Wall -O3 -fopenmp -I$(INCL) -funroll-loops -DNOMAT -DINTEL -msse4.2
	LIBS= -lm
	#OPTS=-DGNU
endif

SRC=.

OBJS = props2D.o

%.o: $(SRC)/%.c 
	$(CC) -c $< $(CFLAGS) $(OPTS) -fpic

all: libprops

libprops: $(OBJS)
	$(CC) -o $(BIN)/$@.so $^ $(CFLAGS) $(LIBS) -shared
	#ar rcs $@.a $^
	rm -f *.o

clean:
	$(RM) *.o $(BIN)/libprops.so

cleandata:
	$(RM) ill_SEGsaltmodel.npy ill_synthmodel2D.npy traces_SEGsaltmodel.npy traces_synthmodel2D.npy waves_SEGsaltmodel.npy waves_synthmodel2D.npy image_SEGsaltmodel.npy image_synthmodel2D.npy rill_SEGsaltmodel.npy rill_synthmodel2D.npy sill_SEGsaltmodel.npy sill_synthmodel2D.npy
