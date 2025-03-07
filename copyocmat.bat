REM Batch Datei

@ECHO OFF
ECHO.
ECHO Start compressing and storing actual version of OCMat
ECHO.
set TMPDIR=%TEMP%\ocmat
FOR /F "tokens=2*" %%A IN ('REG QUERY "HKEY_CURRENT_USER\Software\Microsoft\Windows\CurrentVersion\Explorer\Shell Folders" /v Personal ^| FIND "REG_SZ"')DO (SET MYDOCUMENTDIR=%%B)

REM change default Matlab directory
set MATLABDIR=%MYDOCUMENTDIR%\Matlab
set OCMATDIR=%CD%
set COPYINFOFILE=%TMPDIR%\copyinfo.txt

REM MAKECAB generates a file ~.rpt, in the first line MAKECAB is written together with the actual date. The tokens of this line are used to get MYDATE, etc.
PUSHD "%TEMP%"
MAKECAB /D RptFileName=~.rpt /D InfFileName=~.inf /f NUL > NUL
FOR /F "tokens=3-7" %%a in ('FIND /I "MAKECAB"^<~.rpt') DO (
    SET "MYDATE=%%c-%%b-%%e"
    SET "MYTIME=%%d"
)
DEL ~.*
POPD
ECHO %TEMP%

IF NOT EXIST "%MATLABDIR%" GOTO WARNINGEND
mkdir %TMPDIR%

robocopy  "%OCMATDIR%" %TMPDIR% /e /xf *  > log:nul
REM xcopy /T/E/Y/C "%OCMATDIR%" %TMPDIR%
REM PAUSE
rd /S/Q "%TMPDIR%\model\usermodel"
mkdir "%TMPDIR%\model\usermodel\out"

ECHO.>%COPYINFOFILE%
ECHO Start compressing and storing actual version of OCMat>>%COPYINFOFILE%
ECHO.>>%COPYINFOFILE%
ECHO Files for options>>%COPYINFOFILE%
xcopy /F/Y "%OCMATDIR%\options" %TMPDIR%\options>>%COPYINFOFILE%
ECHO.>>%COPYINFOFILE%
ECHO Documentation files>>%COPYINFOFILE%
xcopy /F/Y "%OCMATDIR%\doc" %TMPDIR%\doc>>%COPYINFOFILE%
ECHO.>>%COPYINFOFILE%
ECHO Continuation files>>%COPYINFOFILE%
xcopy /F/Y/E "%OCMATDIR%\cont" %TMPDIR%\cont>>%COPYINFOFILE%
ECHO.>>%COPYINFOFILE%
ECHO Continuation files>>%COPYINFOFILE%
xcopy /F/Y/E "%OCMATDIR%\calc" %TMPDIR%\calc>>%COPYINFOFILE%
ECHO.>>%COPYINFOFILE%
ECHO General tools>>%COPYINFOFILE%
xcopy /F/Y/E "%OCMATDIR%\tools" %TMPDIR%\tools>>%COPYINFOFILE%
xcopy /F/Y/E "%OCMATDIR%\class" %TMPDIR%\class>>%COPYINFOFILE%
ECHO.>>%COPYINFOFILE%
ECHO Example model>>%COPYINFOFILE%
xcopy /F/Y "%OCMATDIR%\model\example\harvest2D" %TMPDIR%\model\example\harvest2D>>%COPYINFOFILE%
ECHO.>>%COPYINFOFILE%
ECHO Initialization files>>%COPYINFOFILE%
xcopy /F/Y "%OCMATDIR%\model\initfiles" %TMPDIR%\model\initfiles>>%COPYINFOFILE%
ECHO.>>%COPYINFOFILE%
ECHO Deafault file generation files>>%COPYINFOFILE%
xcopy /F/Y "%OCMATDIR%\model\default\filegeneration" %TMPDIR%\model\default\filegeneration>>%COPYINFOFILE%
ECHO Standard model files>>%COPYINFOFILE%
xcopy /F/Y "%OCMATDIR%\model\standardmodel" %TMPDIR%\model\standardmodel>>%COPYINFOFILE%
xcopy /F/Y "%OCMATDIR%\model\standardmodel\filegeneration" %TMPDIR%\model\standardmodel\filegeneration>>%COPYINFOFILE%
xcopy /F/Y "%OCMATDIR%\model\standardmodel\templatefiles" %TMPDIR%\model\standardmodel\templatefiles>>%COPYINFOFILE%
ECHO.>>%COPYINFOFILE%
ECHO Spatial model files>>%COPYINFOFILE%
xcopy /F/Y "%OCMATDIR%\model\ppdemodel" %TMPDIR%\model\ppdemodel>>%COPYINFOFILE%
xcopy /F/Y "%OCMATDIR%\model\ppdemodel\filegeneration" %TMPDIR%\model\ppdemodel\filegeneration>>%COPYINFOFILE%
xcopy /F/Y "%OCMATDIR%\model\ppdemodel\templatefiles" %TMPDIR%\model\ppdemodel\templatefiles>>%COPYINFOFILE%
ECHO.>>%COPYINFOFILE%
ECHO ODE model files>>%COPYINFOFILE%
xcopy /F/Y "%OCMATDIR%\model\odemodel" %TMPDIR%\model\odemodel>>%COPYINFOFILE%
xcopy /F/Y "%OCMATDIR%\model\odemodel\filegeneration" %TMPDIR%\model\odemodel\filegeneration>>%COPYINFOFILE%
xcopy /F/Y "%OCMATDIR%\model\odemodel\templatefiles" %TMPDIR%\model\odemodel\templatefiles>>%COPYINFOFILE%
ECHO.>>%COPYINFOFILE%
ECHO Impulsemodel model files>>%COPYINFOFILE%
xcopy /F/Y "%OCMATDIR%\model\impulsemodel" %TMPDIR%\model\impulsemodel>>%COPYINFOFILE%
xcopy /F/Y "%OCMATDIR%\model\impulsemodel\filegeneration" %TMPDIR%\model\impulsemodel\filegeneration>>%COPYINFOFILE%
xcopy /F/Y "%OCMATDIR%\model\impulsemodel\templatefiles" %TMPDIR%\model\impulsemodel\templatefiles>>%COPYINFOFILE%
ECHO.>>%COPYINFOFILE%
ECHO Difference model files>>%COPYINFOFILE%
xcopy /F/Y "%OCMATDIR%\model\standarddiffmodel" %TMPDIR%\model\standarddiffmodel>>%COPYINFOFILE%
xcopy /F/Y "%OCMATDIR%\model\standarddiffmodel\filegeneration" %TMPDIR%\model\standarddiffmodel\filegeneration>>%COPYINFOFILE%
xcopy /F/Y "%OCMATDIR%\model\standarddiffmodel\templatefiles" %TMPDIR%\model\standarddiffmodel\templatefiles>>%COPYINFOFILE%
ECHO.>>%COPYINFOFILE%
ECHO Open-loop game files>>%COPYINFOFILE%
xcopy /F/Y "%OCMATDIR%\model\standardopenloopgame" %TMPDIR%\model\standardopenloopgame>>%COPYINFOFILE%
xcopy /F/Y "%OCMATDIR%\model\standardopenloopgame\filegeneration" %TMPDIR%\model\standardopenloopgame\filegeneration>>%COPYINFOFILE%
xcopy /F/Y "%OCMATDIR%\model\standardopenloopgame\templatefiles" %TMPDIR%\model\standardopenloopgame\templatefiles>>%COPYINFOFILE%
ECHO.>>%COPYINFOFILE%
ECHO Documentation files>>%COPYINFOFILE%
xcopy /F/Y/E /EXCLUDE:excludelist.txt "%OCMATDIR%\doc" %TMPDIR%\doc>>%COPYINFOFILE%
ECHO.>>%COPYINFOFILE%
ECHO Files in base directory>>%COPYINFOFILE%
xcopy /F/Y "%OCMATDIR%" %TMPDIR%>>%COPYINFOFILE%

REM if input argument is the wildcard * copy all user specific models
if "%1"=="*" GOTO COPYALLMODEL

REM run through all specified model names of the argument list
FOR %%A IN (%*) DO (
ECHO.>>%COPYINFOFILE%
ECHO Specific model: %%A>>%COPYINFOFILE%
mkdir %TMPDIR%\model\usermodel\%%A
xcopy /F/Y/E /EXCLUDE:excludelist.txt "%OCMATDIR%\model\usermodel\%%A" %TMPDIR%\model\usermodel\%%A>>%COPYINFOFILE%
)
GOTO ZIPNMOVE

:COPYALLMODEL 
ECHO.>>%COPYINFOFILE%
ECHO All specific model files>>%COPYINFOFILE%
xcopy /F/Y/E "%OCMATDIR%\model\usermodel" %TMPDIR%\model\usermodel>>%COPYINFOFILE%
ECHO.>>%COPYINFOFILE%
ECHO All documentation files>>%COPYINFOFILE%
xcopy /F/Y/E  "%OCMATDIR%\doc" %TMPDIR%\doc>>%COPYINFOFILE%
)

:ZIPNMOVE
REM zip into matlab folder
IF "%PROCESSOR_ARCHITEW6432%" == "AMD64" GOTO ZIPNMOVE64BIT
"%PROGRAMFILES%\7-Zip\7z" a -tzip "%MATLABDIR%\old\ocmat%MYDATE%.zip" %TMPDIR% -mx9 > NUL

GOTO CLEAN

:ZIPNMOVE64BIT
"%ProgramW6432%\7-Zip\7z" a -tzip "%MATLABDIR%\old\ocmat%MYDATE%.zip" %TMPDIR% -mx9 > NUL

:CLEAN
rd /S/Q %TMPDIR%
ECHO.
ECHO OCMat is stored in:
ECHO	 "%MATLABDIR%\old\ocmat%MYDATE%.zip"
GOTO END

:WARNINGEND
ECHO Matlab directory "%MATLABDIR%" does not exist. 
ECHO Copying aborted.
ECHO.
:END