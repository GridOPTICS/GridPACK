FC = gfortran
CC = gcc

fc_test.x : fc_test.o overwrite.o
	$(FC) -g fc_test.o overwrite.o -o fc_test.x

fc_test.o : fc_test.F90
	$(FC) -g fc_test.F90 -c

overwrite.o : overwrite.c
	$(CC) -g overwrite.c -c

fc2_test.x : fc2_test.o overwrite2.o
	$(FC) -g fc2_test.o overwrite2.o -o fc2_test.x

fc2_test.o : fc2_test.F90
	$(FC) -g fc2_test.F90 -c

overwrite2.o : overwrite2.c
	$(CC) -g overwrite2.c -c

f_from_c.x : c_main.o for_print.o
	$(FC) -g c_main.o for_print.o -o f_from_c.x

c_main.o : c_main.c
	$(CC) -g c_main.c -c

for_print.o : for_print.F90
	$(FC) -g for_print.F90 -c
