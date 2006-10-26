cd src
for %%i in (*.c) do wpp386 -i=..\include %%i
cd ..\lib
for %%i in (..\src\*.obj) do wlib ..\lib\nl.watcom.lib -+%%i
cd ..\examples
for %%i in (*.c) do wpp386 -i=..\include %%i
for %%i in (*.obj) do wlink file %%i, ..\lib\nl.watcom.lib

cd..

del src\*.obj
del examples\*.obj
