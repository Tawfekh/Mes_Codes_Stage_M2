
nom_inpout=$1

cat ${nom_inpout} | cut -d' ' -f2- | sed -e 's/ /,/g' | sed '1d'| sed -e 's/-/N/g' > $(basename -s .txt ${nom_inpout})_h12_1.txt

echo "done"
