path="../samples"

# grep "Direct configuration" XDATCAR | wc -l
for i in {1..4}
do
index="0000$i"
dir="sample.${index}"

truncate -s 0 "XDATCAR-${index}" 
# This command truncates the file  to a size of 0, 
# effectively removing all content from the file 
# or creating an empty file if it doesn't exist

vasp_xdacar="XDATCAR"

n=0 # numbe of time steps

for j in {1..4}
do
file="${vasp_xdacar}-${j}"
#first file
if [ "$j" -eq 1 ]
then

echo "Processing file: ${file}"
cat ${path}/${dir}/${file} > XDATCAR-${index}

m=`grep "Direct configuration" ${path}/${dir}/${file} | wc -l`
n=`echo "$n +$m" | bc -l`

else
#other files
echo "Processing file: ${file}"
tail -n +8 ${path}/${dir}/${file} > temp.dat

m=`grep "Direct configuration" ${path}/${dir}/${file} | wc -l`
n=`echo "$n +$m" | bc -l`

cat temp.dat >> XDATCAR-${index}
rm temp.dat

fi
done

# Last files
if test -f ${path}/${dir}/${vasp_xdacar} 
then
echo "Processing file: ${file}"
tail -n +8 ${path}/${dir}/${vasp_xdacar} > temp.dat

m=`grep "Direct configuration" ${path}/${dir}/${vasp_xdacar} | wc -l`
n=`echo "$n +$m" | bc -l`

cat temp.dat >> XDATCAR-${index} 
rm temp.dat
fi

m=`grep "Direct configuration" XDATCAR-${index} | wc -l`

if [ "$n" -ne "$m" ];then
echo " Files are not properly merge"
echo " ****Please modify the codes****"
echo " Code will STOP here"
exit
else
echo " Total numbe Time Step run = ${n} for job ${i}"
fi
done

fout="XDATCAR"
truncate -s 0 ${fout}

n=0
for ((j = 1; j <= i; j++)); do

fin="XDATCAR-0000${j}"
#outer-loop
if test -f ${fin}; then
# inner-loop
if [ "${j}" -eq 1 ]; then
m=`grep "Direct configuration" ${fin} | wc -l`
n=`echo "$n +$m" | bc -l`
echo " Total numbe Time Step run = ${m}"
cat ${fin} >>${fout}

else 
m=`grep "Direct configuration" ${fin} | wc -l`
n=`echo "$n +$m" | bc -l`
echo " Total numbe Time Step run = ${m}"
tail -n +8 ${fin} > temp.dat 
cat temp.dat >>${fout}
rm temp.dat
fi # end inner-loop

else
echo "*****${fin} does not exit****** "
fi

done

m=`grep "Direct configuration" ${fout} | wc -l`

if [ "$n" -ne "$m" ];then
echo " Files are not properly merge to ${fout}"
echo " ****Please modify the codes****"
else
echo " Total numbe Time Step run = ${n} for job ${i}"
fi

