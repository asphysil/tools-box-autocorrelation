for i in 1000 2000 4000 6000 8000 
do
cat >input-setup.dat<<!
12022 # number of steps
800   # No. of data skiped
$i  # max=((number of steps) - (number of data skip))/2
2   #Time steps
2   # Type of atoms
4 8  # number of atoms for each type 
178.49 15.999 # atomic mass
XDATCAR
born-charge.dat
!
./exe_correl.x
mv md-ir-spectra.dat md-ir-spectra-$i.dat
mv md-ph-dos.dat md-ph-dos-$i.dat
done
