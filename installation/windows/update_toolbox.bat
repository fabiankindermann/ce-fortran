:: A SHELL SCRIPT FOR UPDATING THE TOOLBOX ON A WINDOWS SYSTEM
::
:: ATTENTION: Fortran must already be installed using our original installation files.
::
:: This code is published under the GNU General Public License v3
::                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
::
:: Author: Fabian Kindermann (contact@ce-fortran.com)

@ECHO off
chcp 1252


SET "location=C:\cygwin"


:: GET INTO ADMINSTRATOR MODE
:checkPrivileges 
NET FILE 1>NUL 2>NUL
if '%errorlevel%' == '0' ( goto gotPrivileges ) else ( goto getPrivileges ) 

:getPrivileges 
if '%1'=='ELEV' (shift & goto gotPrivileges)  

setlocal DisableDelayedExpansion
set "batchPath=%~0"
setlocal EnableDelayedExpansion
ECHO Set UAC = CreateObject^("Shell.Application"^) > "%temp%\OEgetPrivileges.vbs" 
ECHO UAC.ShellExecute "!batchPath!", "ELEV", "", "runas", 1 >> "%temp%\OEgetPrivileges.vbs" 
%SystemRoot%\System32\CScript.exe "%temp%\OEgetPrivileges.vbs" 
exit /B 

:gotPrivileges
@setlocal enableextensions
@cd /d "%~dp0"


:: ASK FOR INSTALLATION DIRECTORY CHANGES
ECHO.
ECHO This script updates our toolbox into the directory: 
ECHO.
ECHO    %location%
ECHO.
ECHO To change the installation directory, read the online instructions.
ECHO.
ECHO ATTENTION: Fortran must already be installed using our original installation files.
ECHO. 
set /p temp=Do you want to continue (y/n)?
ECHO.

set exitres=T
if %temp:~-0,1%==y set exitres=F
if %temp:~-0,1%==Y set exitres=F
if %exitres%==T exit


:: INSTALL THE TOOLBOX
gfortran -c -Werror -fopenmp -Wno-unused -ffree-line-length-none -fimplicit-none -Wall -fcheck=bound,do -ffpe-trap=invalid,zero,overflow -frecursive -g ./../toolbox/toolbox.f90 -o toolbox_debug.o
gfortran -c -O3 -fopenmp -ffree-line-length-none ./../toolbox/toolbox.f90 -o toolbox.o
mkdir "%location%\include\" 2>nul
del /Q "%location%\include\toolbox.mod" 2>nul
del /Q "%location%\include\toolbox.o" 2>nul
del /Q "%location%\include\toolbox_debug.o" 2>nul
move "toolbox.mod" "%location%\include\"
move "toolbox.o" "%location%\include%\"
move "toolbox_debug.o" "%location%\include%\"


:: IF EVERYTHING RAN CORRECTLY, AT THIS POINT EVERYTHING SHOULD BE INSTALLED PROPERLY
:theend
ECHO. 
ECHO. 
ECHO. 
ECHO ...TOOLBOX UPDATE COMPLETED.
ECHO.
ECHO.
ECHO In case you encountered any problem, check on www.ce-fortran.com for help.
ECHO.
ECHO.
set /p temp=Press RETURN to end...
ECHO.
