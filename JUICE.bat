@echo off
rem JUICE launcher: opens a NEW, independent window each time it is run (double-click or
rem from a terminal) and returns immediately, so several experiments can be worked on
rem in parallel.  Uses the venv's pythonw.exe (no console window) if present.
setlocal
set "HERE=%~dp0"
set "PYW=%HERE%.venv\Scripts\pythonw.exe"
if not exist "%PYW%" set "PYW=%HERE%..\..\.venv\Scripts\pythonw.exe"
if not exist "%PYW%" set "PYW=pythonw.exe"
start "" "%PYW%" "%HERE%juice_gui.py" %*
endlocal
