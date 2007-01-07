rem This does noty work
md tmp
icl /Iinclude /Fotmp\ /c src\*.c
lib /OUT:lib\nl.lib tmp\*.obj /LIBPATH:"C:\Program Files\Intel\CPP\Compiler80\Ia32\Lib\"
for %%i in (examples\*.c) do icl /O2 /Iinclude /Fotmp\ /Feexamples\ %%i lib\nl.lib  libguide40.lib mkl_c.lib /link /LIBPATH:"C:\Program Files\Intel\CPP\Compiler80\Ia32\Lib\"
del tmp\*.obj
