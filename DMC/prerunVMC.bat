@echo off

cd ..\VMC
gcc -o singlerunVMC.exe singlerunVMC.c 

@REM gcc -o singlerun.exe singlerun.c

set Nt=%1
set Nw=%2
set Nb=%3
set NbSkip=%4

REM Run the compiled C program
singlerunVMC.exe %Nt% %Nw% %Nb% %NbSkip%




