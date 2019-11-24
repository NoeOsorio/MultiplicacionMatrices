PROGNAME := check_memory

main: memory.h
	mpicc $(PROGNAME).c memory.c -o $(PROGNAME);

clean:
	rm -f $(PROGNAME) *.o
