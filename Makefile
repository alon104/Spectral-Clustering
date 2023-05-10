CFlags = -ansi -Wall -Wextra -Werror -pedantic-errors -lm
spkmeans: spkmeans.o spkmeans.h 
	gcc -o spkmeans spkmeans.o $(CFlags)

spkmeans.o: spkmeans.c 
	gcc -c spkmeans.c $(CFlags)

 