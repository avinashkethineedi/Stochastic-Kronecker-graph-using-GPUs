CFLAGS = -O3 -lm
CC = gcc
MCC = mpic++

all: build

build: skg

skg: skg.o
	$(MCC) $(CFLAGS) $< -o $@

skg.o: skg.cpp
	$(MCC) $(CFLAGS) $< -c -o $@

clean:
	rm -rf skg core* *.o