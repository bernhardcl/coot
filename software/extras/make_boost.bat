set PATH=C:\MinGW\bin;%PATH%
set PREFIX=%install_top_dir%
cd tools\build
pwd
bootstrap.bat gcc
.\b2 install --prefix=%PREFIX% toolset=gcc
set PATH=%PREFIX%\bin;%PATH%
cd ..\..\
b2 -j3 install --prefix=%PREFIX% toolset=gcc --layout=system --with-python --with-regex --with-system --with-serialization link=shared release
REM this is system path, may be problematic
REM b2 -j3 install --prefix=%PREFIX% toolset=gcc --layout=system --with-python --with-regex --with-system --with-serialization link=shared release
REM without debug
REM b2 -j3 install --prefix=%PREFIX% toolset=gcc --with-python --with-regex --with-system --with-serialization link=shared release
