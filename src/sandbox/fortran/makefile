FC = gfortran
CC = gcc

fc_test.x : f_test.o overwrite.o
	$(FC) -g f_test.o overwrite.o -o fc_test.x

f_test.o : f_test.F90
	$(FC) -g f_test.F90 -c

overwrite.o : overwrite.c
	$(CC) -g overwrite.c -c

fc2_test.x : f2_test.o overwrite2.o
	$(FC) -g f2_test.o overwrite2.o -o fc2_test.x

f2_test.o : f2_test.F90
	$(FC) -g f2_test.F90 -c

overwrite2.o : overwrite2.c
	$(CC) -g overwrite2.c -c
