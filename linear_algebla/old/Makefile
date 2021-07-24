# Makefile

main: main.o math_lib.o math_lib.h
		g++ -std=gnu++17 -Wall -Wextra -O2 -o main main.o math_lib.o
		
main.o: main.cpp math_lib.h 
		g++ -std=gnu++17 -Wall -Wextra -O2 -o main.o -c main.cpp

math_lib.o: math_lib.cpp math_lib.h
		g++ -std=gnu++17 -Wall -Wextra -O2 -o math_lib.o -c math_lib.cpp




clean:
		rm -f main *.o
