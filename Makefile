CC=gcc
CFLAGS=-lm

mapping_CC: mapping_CC.cpp
	$(CC) -o mapping_CC mapping_CC.cpp $(CFLAGS) 
