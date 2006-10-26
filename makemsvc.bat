md tmp
cl /O2 /Iinclude /Fotmp\ /c src\*.c
lib /OUT:lib\nl.lib tmp\*.obj
for %%i in (examples\*.c) do cl /O2 /Iinclude /Fotmp\ /Feexamples\ %%i lib\nl.lib
del tmp\*.obj
rmdir tmp

