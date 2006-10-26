for %%i in (*.out) do call lgrind.bat -i -o %%i.tex  %%i
for %%i in (*.c) do call lgrind.bat -i -o %%i.tex  %%i
