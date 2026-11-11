@echo off
rem JUICE launcher: opens the JUICE start menu (new experiment, saved experiment, new study,
rem results of one experiment, results of a study).  Every run is a NEW, independent menu and
rem returns immediately, so another window can be opened while others are still open.
rem A file dropped on this launcher (or its shortcut) opens straight in the right window:
rem   .cs -> results viewer, .study.json -> study window, .experiment.json -> experiment window.
rem Uses the venv's pythonw.exe (no console window) if present.
setlocal
set "HERE=%~dp0"
set "PYW=%HERE%.venv\Scripts\pythonw.exe"
if not exist "%PYW%" set "PYW=%HERE%..\..\.venv\Scripts\pythonw.exe"
if not exist "%PYW%" set "PYW=pythonw.exe"
rem the 'juice' package is in this folder
set "PYTHONPATH=%HERE:~0,-1%;%PYTHONPATH%"
start "" "%PYW%" -m juice.gui.juice_start %*
endlocal
