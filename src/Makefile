all:
	gcc -lm -g -Wall -std=c99 -c md.c
	gcc -lm -g -Wall -std=c99 -c verlet.c
	gcc -lm -g -Wall -std=c99 -c cell_list.c
	gcc -lm -g -Wall -std=c99 -c lj.c
	gcc md.o verlet.o cell_list.o lj.o -o md
