@set LANG=en
title Coot's findwaters


set COOT_PREFIX=%~dp0\..

set COOT_HOME=%USERPROFILE%\COOT
set COOT_BACKUP_DIR=%COOT_HOME%\coot-backup

set COOT_SHARE=%COOT_PREFIX%\share

REM shouldnt need the dictionary
REM if not exist "%CLIBD_MON%" (
REM   echo no $CLIBD_MON found setting COOT_REFMAC_LIB_DIR
REM   set COOT_REFMAC_LIB_DIR=%COOT_SHARE%\coot\lib
REM )

set COOT_SCHEME_DIR=%COOT_SHARE%/coot/scheme
set COOT_STANDARD_RESIDUES=%COOT_SHARE%\coot\standard-residues.pdb
set COOT_PIXMAPS_DIR=%COOT_SHARE%\coot\pixmaps
set COOT_RESOURCES_FILE=%COOT_SHARE%\coot\cootrc
set COOT_DATA_DIR=%COOT_SHARE%\coot
set COOT_REF_STRUCTS=%COOT_SHARE%\coot\reference-structures
set COOT_PYTHON_DIR=%COOT_PREFIX%\lib\python3.11\site-packages\coot
set PYTHONHOME=%COOT_PREFIX%

set SYMINFO=%COOT_SHARE%\coot\syminfo.lib

set PATH=%COOT_PREFIX%\bin;%COOT_PREFIX%\lib;%PATH%

%~n0-bin.exe %*
