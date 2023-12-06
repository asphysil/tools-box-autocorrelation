
FC=gfortran
FCFLAG= -Wall -Wextra -pedantic -fimplicit-none -fbacktrace
#-Wl,-ld_classic -Wextra -pedantic -fimplicit-none -fbacktrace 
LIBPATH=
LIB= ${LIBPATH} 
#-llapack -lblas
objects= correlation_cal.o\
         sys_dtype.o\
         cal_properties.o\
         time_correlation_cal.o\
         readfiles_list.o writefiles_list.o\
         main.o
        

exe_correl.x:$(objects) 
	$(FC) ${FCFLAG} -o exe_correl.x $(objects) $(LIB) 
correlation_cal.o: correlation_cal.f90
	${FC} ${FCFLAG} -c correlation_cal.f90
sys_dtype.o: sys_dtype.f90
	${FC} ${FCFLAG} -c sys_dtype.f90 
cal_properties.o: cal_properties.f90
	${FC} ${FCFLAG} -c cal_properties.f90 
time_correlation_cal.o: time_correlation_cal.f90
	${FC} ${FCFLAG} -c time_correlation_cal.f90
readfiles_list.o: readfiles_list.f90 
	${FC} ${FCFLAG} -c readfiles_list.f90 
writefiles_list.o: writefiles_list.f90 
	${FC} ${FCFLAG} -c writefiles_list.f90 
main.o: main.f90
	$(FC) ${FCFLAG} -c main.f90 $(LIB) 
clean:
	rm $(objects)  *.mod *.o exe_correl.x
