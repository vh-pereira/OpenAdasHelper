@echo off

python --version >nul 2>&1
IF ERRORLEVEL 1 (
    echo Python is not installed. Please install it from https://www.python.org/downloads/
    exit /b 1
)

IF NOT EXIST "venv" (
    echo Creating virtual environment...
    python -m venv .venv
) ELSE (
    echo Virtual environment already exists.
)

call .venv\Scripts\activate

IF EXIST "requirements.txt" (
    echo Installing dependencies...
    pip install -r requirements.txt
    IF ERRORLEVEL 1 (
        echo Error installing dependencies. Deactivating venv...
        call venv\Scripts\deactivate.bat
        exit /b 1
    )
)

