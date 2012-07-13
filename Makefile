CFLAGS=-O3 -Wall
LDFLAGS=-lpthread -g -lm

run.exe:  
	    icc -o run.exe extended-shinkage.c $(LDFLAGS) $(CFLAGS) 

clean: 
		rm run.exe

