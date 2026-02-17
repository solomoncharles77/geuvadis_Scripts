

cd geuvadis_Scripts/

split -l 30 allSeq2.txt bamList_ -d -a 2 && 
for file in bamList_*; do 
  mv "$file" "${file}.txt"
done



