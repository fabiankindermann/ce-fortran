:: A SHELL SCRIPT FOR UNINSTALLING FORTRAN FROM A WINDOWS SYSTEM
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


:: ASK FOR UNINSTALLATION PROCESS

ECHO.
ECHO This script completely uninstalls Chocolatey/Fortran/GNU Plot/Geany from your system location:
ECHO.
ECHO    %location%
ECHO.
set /p temp=Do you want to continue (y/n)?
ECHO.

set exitres=T
if %temp:~-0,1%==y set exitres=F
if %temp:~-0,1%==Y set exitres=F
if %exitres%==T exit


:: UNINSTALL GEANY
choco uninstall -y geany

:: UNINSTALL gnuplot
choco uninstall -y gnuplot

:: UNINSTALL CYG-GET
choco uninstall -y cyg-get

:: UNINSTALL CYGWIN
choco uninstall -y cygwin

:: UNINSTALL CHOCOLATEY
rmdir /S /Q %ChocolateyInstall%

:: UPDATE REGISTRY

:: backup current registry key
reg add "HKEY_LOCAL_MACHINE\SYSTEM\CurrentControlSet\Control\Session Manager\Environment" /v Path_backup /t REG_EXPAND_SZ /d "%PATH%" /f

:: update path variable (if necessary)
SetLocal EnableDelayedExpansion

set "unwanted=;%location%\bin"
ECHO.%PATH% | findstr /C:"%location%\bin" 1>nul
if errorlevel 0 (
    set "PATH=!PATH:%unwanted%=!"
)

set "unwanted=;%location%\usr\bin"
ECHO.%PATH% | findstr /C:"%location%\usr\bin" 1>nul
if errorlevel 0 (
    set "PATH=!PATH:%unwanted%=!"
)


set "unwanted=;%ChocolateyInstall%\bin"
ECHO.%PATH% | findstr /C:"%ChocolateyInstall%\bin" 1>nul
if errorlevel 0 (
    set "PATH=!PATH:%unwanted%=!"
)

reg add "HKEY_LOCAL_MACHINE\SYSTEM\CurrentControlSet\Control\Session Manager\Environment" /v Path /t REG_EXPAND_SZ /d "%PATH%" /f
reg query "HKEY_LOCAL_MACHINE\SYSTEM\CurrentControlSet\Control\Session Manager\Environment"

EndLocal

reg delete "HKEY_LOCAL_MACHINE\SYSTEM\CurrentControlSet\Control\Session Manager\Environment" /v ChocolateyInstall /f


:: IF EVERYTHING RAN CORRECTLY, AT THIS POINT EVERYTHING SHOULD BE UNINSTALLED PROPERLY

ECHO. 
ECHO. 
ECHO. 
ECHO ...UNINSTALLATION PROCESS COMPLETED.
ECHO.
ECHO.
ECHO In case you encountered any problem, check on www.ce-fortran.com for help.
ECHO.
ECHO.
set /p temp=Press RETURN to end...
ECHO.
