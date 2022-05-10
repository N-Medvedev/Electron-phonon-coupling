@echo off

:: read argument from user
    set arg1=%1 
    
::  in case of empty argument, assume no debug
    if "%1"=="" ( 
        set arg1=NODEBUG
        )

:: choose how to compile, for debug or for release
    IF %arg1%==DEBUG (
        ifort.exe /debug:all /check:all /check:bounds /fp:precise /fpe-all:0 /Qopenmp /Qftz- /Qfp-stack-check /Od /Zi /traceback /gen-interfaces /warn:all /warn:nounused /Qvec-report1 /fpp /Qtrapuv /dbglibs Analyzing_DOS_main.f90 -o Analyzing_DOS.exe /link /stack:9999999999 

    ) ELSE (
        ifort.exe /F9999999999 /O3 /Qipo /Qvec-report1 /fpp /Qopenmp /heap-arrays Analyzing_DOS_main.f90 -o Analyzing_DOS.exe /link /stack:9999999999 
    )

:: Remove files that are no longer needed
del *.obj *.mod  
:: *.dat


