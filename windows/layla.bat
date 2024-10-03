@set LANG=en
title LAYLA


set COOT_PREFIX=%~dp0\..

set COOT_SHARE=%COOT_PREFIX%\share
set COOT_DATA_DIR=%COOT_SHARE%\coot

set PATH=%COOT_PREFIX%\bin;%PATH%

%~n0.exe %*
