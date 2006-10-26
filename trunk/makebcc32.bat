cd src
bcc32 -O2 -c -I..\include *.c
cd ..\lib
for %%i in (..\src\*.obj) do tlib nl.bcc32.lib -+%%i
cd ..\examples
for %%i in (*.c) do bcc32 -O2 -I..\include %%i ..\lib\nl.bcc32.lib
cd ..

del src\*.obj
del examples\*.obj
del examples\*.tds
