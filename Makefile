CPPFLAGS ?= -std=c++11 -g

Test : Test.o fibonacci.o numbers.o lstrings.o geometry.o sort.h numbers.h geometry.h
	g++ $(CPPFLAGS) Test.o fibonacci.o numbers.o lstrings.o geometry.o -o Test
		
clean : 
	rm -f Test *.o core*
	
