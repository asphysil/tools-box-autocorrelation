program main_program
    ! pupose to calculated : phonon DOS, Ir spectra 
    
    use nrtype 
    use sys_constants, only : ndim 
    use  filesystem
    use data_struc
    !use test_data 
    use readfiles, only: read_vasp_xdatcar,read_test_data

    use sys_properties_cal, only: cal_velocity_forward_method,cal_polarization

    use time_correl, only : cal_vel_correl, cal_phonon_dos,cal_pol_correl,cal_ir_spectra
    
    use writefiles, only : write_vasp_md,write_vasp_vac,write_vasp_ph_dos,&
    write_test_fftw, write_vasp_md_pol, write_vasp_pol_acf, write_vasp_ir_spectra
    implicit none

! local 
    real(dp), dimension(:), allocatable :: loc_mass 
    real(dp) :: log_result 
    real(dp) :: xn 
    integer :: nvdata2pow
    integer :: npol2pow
    integer :: nskip
    integer :: nirdata 
    integer :: nvv_corr 
    integer :: npp_corr 

    integer :: i, j, k 
    integer :: ierror
    integer :: junk 
    LOGICAl :: fexist
    character(len=15) :: fborn_chg
    character(len=7) :: fxdatcar

    !!
!!!!!!!!!!!!!!!!!
INQUIRE(FILE="input-setup.dat", exist=fexist)
    IF (fexist) THEN
        PRINT*, "input-setup.dat file exist"
        open(unit=12, file='input-setup.dat', action='read')
    ELSE
      PRINT*, "input-setup.dat file does exist in this folder and program will STOP here"
    STOP
ENDIF

read(12, *) ntot  ! number of md step 
print*, ' number of md step', ntot 
read(12,*) nskip !
print*,'**** Number of data skiped ****', nskip

read(12,*) nirdata 
print*,'**** Number of data will be used for phonon DOS and IR spectra cal. is****', nirdata

read(12, *) dt 
print*, ' time step =', dt 
read(12,*) ntype ! number type of atoms 
print*, ' number of type of atoms', ntype 

allocate(ntype_atms(ntype)) !

    read(12,*) (ntype_atms(i), i=1,ntype)

natms=sum(ntype_atms, dim=1) 

print*, 'total number of atoms', natms
allocate(loc_mass(ntype), atms_mass(natms))

do i=1, ntype
print*, i, " th type of atoms=", ntype_atms(i)
enddo 

read(12,*) (loc_mass(i), i=1, ntype)
k = 0 
do i=1, ntype 
    do j=1, ntype_atms(i)
        k = k + 1
        atms_mass(k) = loc_mass(i)
        print*, 'Atomic Mass ', loc_mass(i), 'type = ', i 
    enddo
enddo 
! reading file name
read(12,*) fxdatcar
read(12,*) fborn_chg
close(12) 
deallocate(loc_mass)

! reading born change 
allocate(born_charg(ndim, ndim, natms))

!!!!!!!!!!!!!!!!!
INQUIRE(FILE="born-charge.dat", exist=fexist)
    IF (fexist) THEN
        PRINT*, "born-charge.dat file exist"
        open(unit=12, file='born-charge.dat', action='read')
    ELSE
      PRINT*, "born-charge.dat file does exist in this folder and program will STOP here"
    STOP
ENDIF
 
print*, '**Reading Born Charge...'

!read(12,*)
read(12,*)

do i =1, natms
    read(12,*)
    do j=1,3 
        read(12,*) junk, (born_charg(j,k,i), k=1,3)
        print*,  (born_charg(j,k,i), k=1,3), 'atoms No', i 
    enddo
    print*, '  '
enddo
print*, '**End reading Born Charge**'


print*, '**** please check first,  how many data to be skip******'
print*, ' Please modify number of data to be used for ACF'
!nskip=600 


ndata = ntot-nskip 
if (ntot<nskip) then 
    print*, ' The number of data set must be greater than 100'
    print*, ' program will stop here'
    stop 
else 
    print*, '*** Number of data to be used for the calculation is *** ', ndata  
endif

nvdata = ndata - 1 ! number of data for velocity calculation 

! acf 
! ! Find the closest lower  powers of 2
xn = real(nvdata, kind=dp)! number of data to be used for ACF
log_result = log(xn) / log(2.0)
nvdata2pow = 2**(int(log_result))
print*, '*** Number of data to be used for the calculation is *** ', nvdata2pow
nc = 2*nvdata2pow ! correlated signal should be 2*nstep points 
! https://github.com/elcorto/pwtools/blob/master/doc/source/written/background/phonon_dos.rst
! for fft 
if ( nirdata >= int(nc/2)) then
! xn = real(2000, kind=dp)
!else
! xn = real(nvdata, kind=dp)
print*, '************Full data set of ACF is used for phonon dos calculation********'
xn = real(int(nc/2), kind=dp)
log_result = log(xn) / log(2.0)
nc_fft = 2**(int(log_result))
else
xn = real(nirdata, kind=dp)
log_result = log(xn) / log(2.0)
nc_fft = 2**(int(log_result)+1)
nvv_corr = nirdata 
endif

if (nirdata > nc/2 ) then 
    print*, ' please decrese number of nirdata '
    print*, ' program will stop here'
    stop 
endif 

!log_result = log(xn) / log(2.0)
!nc_fft = 2**(int(log_result)+1)
print*, '********phonon DOS********'
print*, '**Total number data**', ntot,  '**number of velocity data**', nvdata
print*, '**Colsed integer power of 2 for number **', nvdata,  ' **is** ', nvdata2pow
print*, '**number of correlation data **', nc/2
print*, '**Colsed integer power of 2 for number**', nirdata,  '**is**',  nc_fft
print*, '********END phonon DOS********'
! 
!********************For IR spectra *********************

! acf 
xn = real(ndata, kind=dp)! number of data to be used for ACF 
! ! Find the closest lower  powers of 2
log_result = log(xn) / log(2.0)
npol2pow = 2**(int(log_result))
nc_ir = 2*npol2pow

! for fft 
if ( nirdata >= int(nc_ir/2)) then
! xn = real(2000, kind=dp)
!else
! xn = real(nvdata, kind=dp)
print*, '************Full data set of ACF is used for phonon dos calculation********'
xn = real(int(nc_ir/2), kind=dp)
log_result = log(xn) / log(2.0)
nc_ir_fft = 2**(int(log_result))
else
xn = real(nirdata, kind=dp)
log_result = log(xn) / log(2.0)
nc_ir_fft = 2**(int(log_result)+1)
npp_corr = nirdata
endif

if (nirdata > nc_ir/2 ) then 
print*, ' please decrese number of nirdata '
print*, ' program will stop here'
stop 
endif 

print*, '********IR********'
print*, '**Total number data**', ntot,  '**number of data for polarization calculation**', ndata
print*, '**Colsed integer power of 2 for number **', ndata,  ' **is** ', npol2pow
print*, '**number of correlation data **', nc_ir/2
print*, '**Colsed integer power of 2 for number**', nirdata,  '**is**',  nc_ir_fft
! 
print*, '********END IR********'
!
allocate(latt_vec(ndim, ndim))



allocate(pos(ndim, natms, ntot)) ! Total number data set read

allocate(atms_vel(ndim, natms, nvdata))
allocate(vel_cdiff(ndim, natms, nvdata))
allocate(vel_fdiff(ndim, natms, nvdata))
!allocate(vel_bdiff(ndim, natms, ndata))
allocate(ek(nvdata), temp(nvdata))

allocate(cvel(ndim, natms, nvdata))
allocate(vv_corr(nc))
allocate(ph_dos(nc_fft/2) )

! for IR spectra 

allocate(pol(ndim, ndata))
allocate(pp_corr(nc_ir))
allocate(ir_sp(nc_ir_fft/2))

! close the file

! intialized 
do i=1, nvdata
    ek(i)=0.0
    temp(i)=0.0
enddo




read_fileid=15 
open(unit=read_fileid, file='XDATCAR', action='read')
call read_vasp_xdatcar(read_fileid, ntype,  natms, ntot,  pos)
close(read_fileid)

! Velocity
call  cal_velocity_forward_method(natms, nskip, nvdata) 
!
! ! first check how many data to be skip 
 call cal_vel_correl(nc,  nvdata, vv_corr)

! FFT of auto-correlation
call cal_phonon_dos(nvv_corr, vv_corr(1:nvv_corr), nc_fft, ph_dos)
!
call cal_polarization(natms, nskip, ndata)
call cal_pol_correl(nc_ir, ndata, pp_corr)
call cal_ir_spectra(npp_corr, pp_corr(1:npp_corr), nc_ir_fft, ir_sp)

write_fileid1=20
open(unit=write_fileid1, file='md-vasp.dat', action='write')
call write_vasp_md(write_fileid1, nvdata, dt)
close(write_fileid1)

write_fileid2=21
open(unit=write_fileid2, file='md-vac.dat', action='write')
call write_vasp_vac(write_fileid2, nvdata, dt)
close(write_fileid2)

write_fileid3=22
open(unit=write_fileid3, file='md-ph-dos.dat', action='write')
call write_vasp_ph_dos(write_fileid3, nc_fft/2, dt)
close(write_fileid3)

! for ir-spectra 
write_fileid1=20
open(unit=write_fileid1, file='md-vasp-pol.dat', action='write')
call write_vasp_md_pol(write_fileid1, ndata, dt)
close(write_fileid1)

write_fileid2=21
open(unit=write_fileid2, file='md-pol-acf.dat', action='write')
call write_vasp_pol_acf(write_fileid2, ndata, dt)
close(write_fileid2)

write_fileid3=22
open(unit=write_fileid3, file='md-ir-spectra.dat', action='write')
call write_vasp_ir_spectra(write_fileid3, nc_ir_fft/2, dt)
close(write_fileid3)


! ! ! To test subroutine 
!  ntest=6000
! ! for fft 
! xn = real(ntest, kind=dp)
! log_result = log(xn) / log(2.0)
! nc_fft = 2**(int(log_result)+1)
!  allocate(t_test(ntest), y_test(ntest), x_test(nc_fft), w_test(nc_fft/2))
!  read_fileid = 12
!  open(unit=read_fileid, file='expanded_data.txt', action='read')
!  call read_test_data(read_fileid, ntest, t_test, y_test)
! close(12)
! do i =1, nc_fft
!     x_test(i) = 0.0
! enddo
! !x_test(1:ntest) = y_test(1:ntest)
! do i=1, ntest
!     !print*, y_test(i)
!      x_test(i) = y_test(i)
! enddo
! dt_test = t_test(2)-t_test(1)
! print*, dt_test
! call cal_phonon_dos(nc_fft, x_test, w_test)
!  write_fileid4=23
! open(unit=write_fileid4, file='test-ref.dat', action='write')
! call write_test_fftw(write_fileid4, ntest, nc_fft/2, dt_test, w_test)
! close(write_fileid4)

! close(write_fileid1)
! close(write_fileid2)
! close(write_fileid3)
! close(write_fileid4)

end program main_program
