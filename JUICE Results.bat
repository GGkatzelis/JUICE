@echo off
rem JUICE results viewer: opens a result (.cs) directly, without the start menu.
rem (The JUICE icon / JUICE.bat covers this too: "Look at results" in its start menu.)
rem   double-click            -> "one experiment or a study?", then a file picker
rem   drop a .cs on this file -> that result opens at once (also works on the shortcut)
rem Every call is a new, independent window; uses the venv's pythonw.exe when present.
setlocal
set "HERE=%~dp0"
set "PYW=%HERE%.venv\Scripts\pythonw.exe"
if not exist "%PYW%" set "PYW=%HERE%..\..\.venv\Scripts\pythonw.exe"
if not exist "%PYW%" set "PYW=pythonw.exe"
start "" "%PYW%" "%HERE%results_viewer.py" %*
endlocal
