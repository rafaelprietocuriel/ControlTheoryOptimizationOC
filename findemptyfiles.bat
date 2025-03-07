@ECHO OFF
SET topLevel=%CD%
ECHO. 2>%topLevel%\EmptyFiles.txt
FOR /D /R %%D IN (*) DO (
  CD %%D 
  CALL :innerLoop
)
CD %topLevel%
FOR /F "usebackq delims=" %%D IN (`"DIR /AD/B/S | SORT /R"`) DO (ECHO "%%D" >> %topLevel%\EmptyFiles.txt) 
GOTO :break

:innerLoop
  FOR /F "delims=" %%F IN ('DIR/B/A-D/OS') DO IF %%~zF EQU 0 (ECHO "%CD%\%%F" >> %topLevel%\EmptyFiles.txt) ELSE (GOTO :break)

:break