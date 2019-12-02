#!/bin/bash
file='Makefile'
cd src
 cp $file $file".orig"
 echo have file  $file
 sed -e 's?\/\* Enumerations from ?\/\/\* Enumerations from ?' $file >tmp.txt
 mv tmp.txt $file
cd -

