@echo off

set PROJ=dawg
set PROJ_DISTS=dawg-2*
set MAKE=nmake
set CMAKE=cmake
set SVN=svn
svn info | findstr /b URL | perl -pe "s!^URL: (.+)/releng$!$1!" > url.tmp
set /P REPOS=<url.tmp
del url.tmp

echo.
echo Building distributions for %REPOS% ...

set RELENG_DIR="%TEMP%\%PROJ%-releng.%RANDOM%"
mkdir %RELENG_DIR% || exit /B 1

echo Using temp directory %RELENG_DIR% ...
echo.

set DEST_DIR="%CD%"
set SOURCE_DIR="%RELENG_DIR%\source"
set BUILD_DIR="%RELENG_DIR%\build"

%SVN% co -q %REPOS% %SOURCE_DIR% || exit /B 1

mkdir %BUILD_DIR% || exit /B 1
cd %BUILD_DIR% || exit /B 1

call "C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\vcvarsall.bat" x86

%CMAKE% -G "NMake Makefiles" %SOURCE_DIR% -DCMAKE_BUILD_TYPE=Release -DBoost_USE_STATIC_LIBS=yes -DLIBDAWG_USE_STATIC_LIBS=yes
if %ERRORLEVEL% NEQ 0 goto :end
%MAKE%
if %ERRORLEVEL% NEQ 0 goto :end
%MAKE% package
if %ERRORLEVEL% NEQ 0 goto :end
%MAKE% package_source
if %ERRORLEVEL% NEQ 0 goto :end

echo.
echo Copying distribution packages ...

xcopy /Y %PROJ_DISTS% %DEST_DIR%

echo.
echo Cleaning up ...

cd %DEST_DIR%
rd /S /Q %RELENG_DIR%

:end
cd %DEST_DIR%