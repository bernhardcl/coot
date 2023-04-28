@set LANG=en
title LIDIA


set COOT_PREFIX=%~dp0\..

set COOT_SHARE=%COOT_PREFIX%\share
set COOT_DATA_DIR=%COOT_SHARE%\coot

set PATH=%COOT_PREFIX%\bin;%COOT_PREFIX%\python27;%PATH%

%~n0-bin.exe %*
