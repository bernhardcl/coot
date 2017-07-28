make_file='Makefile'
for file in `find ./* -name ${make_file}` 
do
 cp $file $file".orig"
 echo have file  $file
 sed -e 's?-version-info?-no-undefined -version-info?' $file >tmp.txt
 mv tmp.txt $file
done
