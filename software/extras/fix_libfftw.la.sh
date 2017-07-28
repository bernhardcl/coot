make_file='*fftw*.la'
echo $make_file
for file in `find ./* -name '*fftw*.la'` 
do
 cp $file $file".orig"
 echo have file  $file
 sed -e "s?dependency_libs=.*?dependency_libs=''?" $file >tmp.txt
 mv tmp.txt $file
done
