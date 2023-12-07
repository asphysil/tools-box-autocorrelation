MODULE writefiles
    use nrtype
    use sys_constants, only : ndim, fs2sec
    use data_struc, only: temp, ek, vv_corr, ph_dos
    IMPLICIT NONE 
 
 private    
 public ::  write_vasp_md,write_vasp_vac,write_vasp_ph_dos,&
 write_test_fftw, write_vasp_md_pol, write_vasp_pol_acf, write_vasp_ir_spectra
contains 
SUBROUTINE write_vasp_md(fileID, ndata, dt)
    IMPLICIT NONE 
    INTEGER, INTENT(IN) :: fileID, ndata
    real(dp), intent(in) :: dt

  ! local 
    integer :: i, j
    real(dp) :: t 

  do i=1, ndata
    t = dt*(i-1)
    write(fileID, 200) t, ek(i), temp(i)
  enddo

200 format(1X, 3F15.8)
END SUBROUTINE write_vasp_md

  
SUBROUTINE write_vasp_vac(fileID, ndata, dt)
  IMPLICIT NONE 
  INTEGER, INTENT(IN) :: fileID, ndata
  real(dp), intent(in) :: dt
! local 
  integer :: i, j
  real(dp) :: t 

do i=1, ndata
  t = dt*(i-1)
  write(fileID, 200) t, vv_corr(i)
enddo

200 format(1X, 2F15.8)
END SUBROUTINE write_vasp_vac


SUBROUTINE write_vasp_ph_dos(fileID, ndata, dt)
  IMPLICIT NONE 
  INTEGER, INTENT(IN) :: fileID, ndata
  real(dp), intent(in) :: dt
! local 
  integer :: i, j
  real(dp) :: f0, f

! 
! f0 = 1.0/(dt*N) ! N => Total number data used for FFT and dt is time step
f0 = 1.0/(2*dt* real(ndata, kind=dp))

do i=1, ndata
  f = (i-1)*f0*1.0E3 ! in THz   
  write(fileID, 200) f, ph_dos(i)
enddo

200 format(1X, 2F15.8)
END SUBROUTINE write_vasp_ph_dos


SUBROUTINE write_test_fftw(fileID, norg, ndata, dt_test, y)
  IMPLICIT NONE 
  INTEGER, INTENT(IN) :: fileID, norg, ndata
  real(dp), intent(in) :: dt_test, y(ndata)

! local 
  integer :: i, j
  real(dp) :: f0, f

! 
!print*, 'ok'
! https://aiichironakano.github.io/phys516/VAC.pdf
! f0 = 1.0/(dt*N) ! N => Total number data used for FFT and dt is time step
 !  
f0 = 1.0/(2*dt_test* real(ndata, kind=dp))

do i=1, norg/2 
  f = (i-1)*f0*2.0*PI_D ! omega 
  write(fileID, 300) f, y(i)
enddo

300 format(1X, 2F15.8)
END SUBROUTINE write_test_fftw

subroutine write_vasp_md_pol(fileID, ndata, dt )
  use data_struc, only: pol 
  IMPLICIT NONE 
  INTEGER, INTENT(IN) :: fileID, ndata
  real(dp), intent(in) :: dt
   ! local 
  integer :: i, j
  real(dp) :: t 

do i=1, ndata
  t = dt*(i-1)
  write(fileID, 200) t, (pol(j,i), j=1,3)
enddo
200 format(1X, 4F15.8)
end subroutine write_vasp_md_pol

subroutine write_vasp_pol_acf(fileID, ndata, dt)
  use data_struc, only: pp_corr
  IMPLICIT NONE 
  INTEGER, INTENT(IN) :: fileID, ndata
  real(dp), intent(in) :: dt
   ! local 
  integer :: i, j
  real(dp) :: t 

  do i=1, ndata
    t = dt*(i-1)
    write(fileID, 200) t, pp_corr(i)
  enddo
  
  200 format(1X, 2F15.8)
end subroutine write_vasp_pol_acf

subroutine write_vasp_ir_spectra(fileID, ndata, dt)
  use data_struc, only : ir_sp 
  IMPLICIT NONE 
  INTEGER, INTENT(IN) :: fileID, ndata
  real(dp), intent(in) :: dt
! local 
  integer :: i, j
  real(dp) :: f0, f

! 
! f0 = 1.0/(dt*N) ! N => Total number data used for FFT and dt is time step
f0 = 1.0/(2*dt* real(ndata, kind=dp))

do i=1, ndata
  f = (i-1)*f0*1.0E3 ! in THz   
  write(fileID, 200) f, ir_sp(i)
enddo

200 format(1X, F15.10, 1X, F15.10)
  
end subroutine write_vasp_ir_spectra
END MODULE writefiles
