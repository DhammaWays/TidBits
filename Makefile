CPPFLAGS ?= -std=c++11 -g

Test : Test.o fibonacci.o numbers.o lstrings.o sort.h numbers.h
	g++ $(CPPFLAGS) Test.o fibonacci.o numbers.o lstrings.o -o Test
		
clean : 
	rm -f Test *.o core*
	
