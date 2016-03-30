#FIXME fix help.cc (not al need replacment
file_ls="cprimitive.h
help.cc
help.h
help_globals.h
splineset.cc
texture.h
texture_globals.h"

for file in $file_ls
do
  echo patching $file
  sed -e 's/#if defined (_WIN32)/#if defined (_WIN32) \&\& not defined (WINDOWS_MINGW)/' $file >tmp.txt
  mv tmp.txt $file
done
