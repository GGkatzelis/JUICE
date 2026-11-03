@echo off
rem JUICE results viewer: opens a result (.cs) directly, without the step window.
rem   double-click            -> file picker starting in the runs folder of the TOML
rem   drop a .cs on this file -> that result opens at once (also works on the shortcut)
rem Every call is a new, independent window; uses the venv's pythonw.exe when present.
setlocal
set "HERE=%~dp0"
set "PYW=%HERE%.venv\Scripts\pythonw.exe"
if not exist "%PYW%" set "PYW=%HERE%..\..\.venv\Scripts\pythonw.exe"
if not exist "%PYW%" set "PYW=pythonw.exe"
start "" "%PYW%" "%HERE%results_viewer.py" %*
endlocal
