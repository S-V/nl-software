md tmp
cd tmp
gcc -O3 -Wall -I ../include -c ../src/*.c
ar -ru ../lib/nl.a *.o
for %%i in (..\examples\*.c) do gcc -O3 -Wall -o %%i.exe -I ../include %%i ../lib/nl.a
del *.o
cd ..
