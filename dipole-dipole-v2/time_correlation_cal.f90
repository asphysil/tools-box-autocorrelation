module auto_time_correl_cal 
    use nrtype 
    use sys_constants, only: ndim, Navogadro, ev2J, kb   
    implicit none
    private 
    public :: cal_vv_correl, cal_pp_correl,cal_phonon_dos,cal_ir_spectra
contains

subroutine cal_vv_correl(natms, ndata, nc, vv_correl)
    use data_struc, only :  cvel
    implicit none 
    integer, intent(in) :: natms, ndata, nc
    real(dp), intent(out) :: vv_correl(ndim,nc)

! local 
integer :: n1, n2, tau
integer :: i, j, k, i1, i2
real(dp) :: temp 

if ( (ndata/2) == nc) then
        print*, " Error in 'subroutine cal_velocity_correl' input ndata must be divided by 2"
        print*, "program will stop here"
endif

! intialize
do k=1, 3
vv_correl(k,1:nc)=0.0
enddo

!
! Time correlation
do k= 1, 3
    do j=1, natms

         tau = 0
         do i1 = 1, nc
            n2 = ndata - tau ! max(n2) = ndata/2
            temp = 0.0
            do i2 = 1, n2
               temp = temp + cvel(k, j, i2)*cvel(k, j, i2+tau)
            enddo
            vv_correl(k,i1) = vv_correl(k, i1) + temp
            tau = tau + 1  ! 0*dt, 1*dt, 2*dt ... n*dt 
          enddo

      enddo
    vv_correl(k,1:nc) = vv_correl(k,1:nc)/vv_correl(k,1) ! Normilized 
enddo
end subroutine cal_vv_correl

subroutine cal_pp_correl(ndata, nc, pp_correl)
    use data_struc, only : pol, born_charg 
implicit none 
integer, intent(in) :: ndata, nc
real(dp), intent(out) :: pp_correl(ndim,nc)
!local
integer :: n1, n2, tau 
integer :: i1, i2, k

real(dp) :: temp

if ( (ndata/2) == nc) then
        print*, " Error in 'subroutine cal_pol_correl' input ndata must be divided by 2"
        print*, "program will stop here"
endif
! intialize
do k=1, 3
pp_correl(k,1:nc)=0.0
enddo

! Time correlation
do k= 1, 3
         tau = 0
         do i1 = 1, nc
            n2 = ndata - tau ! max(n2) = nvel/2
            temp = 0.0
            do i2 = 1, n2
               temp = temp + pol(k, i2)*pol(k, i2+tau)
            enddo
            pp_correl(k,i1) = pp_correl(k, i1) + temp
            tau = tau + 1  ! (i1-1)*dt
          enddo
    pp_correl(k,1:nc) = pp_correl(k,1:nc)/pp_correl(k,1) ! Normilized 
enddo

 end subroutine cal_pp_correl 

!!!!!!!!!! Fourier tranform of time correlation function
subroutine cal_phonon_dos(nc, vv_corr, ph_dos)
    USE nrutil, ONLY :assert,assert_eq
    use nr, only :  correl_fft_dp
    IMPLICIT NONE
    integer, intent(in) :: nc 
    REAL(DP), DIMENSION(:), INTENT(IN) :: vv_corr(nc) 
    real(dp), dimension(:), intent(out) :: ph_dos(nc/2)

    COMPLEX(DPC), DIMENSION(size(vv_corr)/2) :: cdata 
    

integer :: i 

  cdata = correl_fft_dp(vv_corr) ! FFT 
 
do i =1, nc/2 
    ph_dos(i) = abs(cdata(i))**2
enddo 

  end subroutine cal_phonon_dos

subroutine cal_ir_spectra(nc, pp_corr, ir_sp)
        USE nrutil, ONLY :assert,assert_eq
        use nr, only :  correl_fft_dp
      IMPLICIT NONE
      integer, intent(in) :: nc 
      REAL(DP), DIMENSION(:), INTENT(IN) :: pp_corr(nc) 
      real(dp), dimension(:), intent(out) :: ir_sp(nc/2)
  
      COMPLEX(DPC), DIMENSION(size(pp_corr)/2) :: cdata 
      
  
  integer :: i 
  
    cdata = correl_fft_dp(pp_corr) ! FFT 
   
  do i =1, nc/2 
      ir_sp(i) = abs(cdata(i))**2
  enddo 
  
    end subroutine cal_ir_spectra   

end module auto_time_correl_cal
