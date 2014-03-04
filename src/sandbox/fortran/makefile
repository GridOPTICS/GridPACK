fc_test.x : f_test.o overwrite.o
	ifort -g f_test.o overwrite.o -o fc_test.x

f_test.o : f_test.F90
	ifort -g f_test.F90 -c

overwrite.o : overwrite.c
	icc -g overwrite.c -c

fc2_test.x : f2_test.o overwrite2.o
	ifort -g f2_test.o overwrite2.o -o fc2_test.x

f2_test.o : f2_test.F90
	ifort -g f2_test.F90 -c

overwrite2.o : overwrite2.c
	icc -g overwrite2.c -c
