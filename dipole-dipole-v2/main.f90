program main_program
    ! pupose to calculated : phonon DOS, Ir spectra 
    
    use nrtype 
    use sys_constants, only : ndim 
    use  filesystem
    use data_struc
    !use test_data 
    use readfiles, only: read_vasp_xdatcar,read_test_data

    use sys_properties_cal, only: cal_velocity_forward_method,cal_polarization

    use auto_time_correl_cal, only : cal_vv_correl, cal_pp_correl
    
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


ndata  = ntot-nskip 
if (ntot<nskip) then 
    print*, ' The number of data set must be greater than 100'
    print*, ' program will stop here'
    stop 
else 
    print*, '*** Number of data to be used for the calculation is *** ', ndata  
endif

nvdata = ndata - 1 ! number of data for velocity calculation 

! acf 
! ! Find the closest number which is divided by 2
! Calculate the closest even number
if (mod(nvdata, 2) == 0) then
    nv_even_no = nvdata
else
    nv_even_no = nvdata  - mod(nvdata, 2)
endif

nc = nv_even_no/2

if ( nirdata >= nc) then
! xn = real(2000, kind=dp)
!else
! xn = real(nvdata, kind=dp)
print*, '************Full data set of ACF is used for phonon dos calculation********'
xn = real(nc, kind=dp)
log_result = log(xn) / log(2.0)
nc_fft = 2**(int(log_result)+1)
else
xn = real(nirdata, kind=dp)
log_result = log(xn) / log(2.0)
nc_fft = 2**(int(log_result)+1)

endif

!log_result = log(xn) / log(2.0)
!nc_fft = 2**(int(log_result)+1)

print*, '**Total number data**', ntot,  '**number of velocity data**', nvdata
print*, '**Colsed integer power of 2 for number **', 2*nvdata,  ' **is** ', nc
print*, '**number of correlation data **', nc
print*, '**Colsed integer power of 2 for number**', nirdata,  '**is**',  nc_fft
! 
!********************For IR spectra *********************

! acf 
! Calculate the closest even number
if (mod(ndata, 2) == 0) then
    nir_even_no = ndata
else
    nir_even_no = ndata  - mod(ndata, 2)
endif

nc_ir = nir_even_no /2 

! for fft 
if ( nirdata >= nc_ir) then
! xn = real(2000, kind=dp)
!else
! xn = real(nvdata, kind=dp)
print*, '************Full data set of ACF is used for phonon dos calculation********'
xn = real(nc_ir, kind=dp)
log_result = log(xn) / log(2.0)
nc_ir_fft = 2**(int(log_result)+1)
else
xn = real(nirdata, kind=dp)
log_result = log(xn) / log(2.0)
nc_ir_fft = 2**(int(log_result)+1)
endif

print*, '**Total number data**', ntot,  '**number of data for polarization calculation**', ndata
print*, '**ndata=', ndata_org,  ' **ndata/2** ', nc_ir
print*, '**number of correlation data **', nc_ir
print*, '**Colsed integer power of 2 for number**', nirdata,  '**is**',  nc_ir_fft
! 

!
allocate(latt_vec(ndim, ndim))



allocate(pos(ndim, natms, ntot)) ! Total number data set read

allocate(atms_vel(ndim, natms, nvdata))
allocate(vel_cdiff(nvdata, natms, ndim))
allocate(vel_fdiff(nvdata, natms, ndim))
!allocate(vel_bdiff(ndim, natms, ndata))
allocate(ek(nvdata), temp(nvdata))

allocate(cvel(ndim, natms, nv_even_no))
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
call cal_velocity_correl(natms, nv_even_no, nc, vv_correl)

! FFT of auto-correlation
call cal_phonon_dos(nc_fft, vv_corr(1:nc_fft), ph_dos)
!
call cal_polarization(natms, nskip, ndata)
call cal_pol_correl(nir_even_no, nc, pp_correl)
call cal_ir_spectra(nc_ir_fft, pp_corr(1:nc_ir_fft), ir_sp)

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
