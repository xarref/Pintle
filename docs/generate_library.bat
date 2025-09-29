@echo off
echo Fluid Property Library Generator
echo ================================
echo.
echo This script generates a comprehensive fluid property library
echo using MATLAB and CoolProp for reference-quality accuracy.
echo.
echo Requirements:
echo   - MATLAB with CoolProp toolbox installed
echo   - CoolProp for MATLAB: https://coolprop.sourceforge.net/coolprop/wrappers/MATLAB/index.html
echo.

pause

echo.
echo Starting MATLAB to generate fluid library...
matlab -batch "generate_fluid_library"

echo.
echo Done! Check for fluid_library.json file (~200KB).
echo Upload this file with your web calculator for reference-quality accuracy.
pause
