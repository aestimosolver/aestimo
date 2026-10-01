@echo off
setlocal
title Aestimo GUI Launcher

echo ==========================================
echo    Aestimo 1D Semiconductor Solver GUI
echo ==========================================
echo.

:: Check if Python is installed
python --version >nul 2>&1
if %errorlevel% neq 0 (
    echo [ERROR] Python not found. Please install Python 3.
    pause
    exit /b 1
)

:: Check if required packages are installed (optional but helpful)
echo [INFO] Verifying dependencies (customtkinter, numpy, matplotlib)...
python -c "import customtkinter, numpy, matplotlib" >nul 2>&1
if %errorlevel% neq 0 (
    echo [WARNING] Some dependencies might be missing.
    echo [INFO] Attempting to install required packages...
    python -m pip install customtkinter numpy matplotlib
)

echo [INFO] Launching Aestimo GUI...
python aestimo_gui.py

if %errorlevel% neq 0 (
    echo.
    echo [ERROR] Aestimo GUI exited with an error.
    pause
)

endlocal
