md tmp
cl /O2 /Iinclude /Fotmp\ /c src\*.c
lib /OUT:lib\nl.lib tmp\*.obj
rem for %%i in (examples\*.c) do cl /O2 /Iinclude /Fotmp\ /Feexamples\ %%i lib\nl.lib libguide40.lib mkl_c.lib
for %%i in (examples\*.c) do cl /O2 /Iinclude /Fotmp\ /Feexamples\ %%i lib\nl.lib libguide40.lib mkl_c.lib
del tmp\*.obj
rmdir tmp

