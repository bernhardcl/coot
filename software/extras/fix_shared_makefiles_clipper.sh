for file in `find ./* -name 'Makefile'` 
do
 cp $file $file".orig"
 sed -e 's?@VERSION_INFO@?-no-undefined @VERSION_INFO@?' $file >tmp.txt
 mv tmp.txt $file
done
