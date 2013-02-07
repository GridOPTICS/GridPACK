See talk.pdf for some information about this project.

convert.F90:  This program converts *MPS files into a binary format that we prefer.

These files support convert.F90:
dict_mod.f90
get_string.f90

These non-PNNL, but open sourced, files are used by convert.F90 to process the dictionary of terms in the MPS file:
dictionary.f90
linkedlist.f90

ftr_andsu7.other.F: This is the actual code that does the linear solving.
Current version alternates "DIIS" and "normal" method.

