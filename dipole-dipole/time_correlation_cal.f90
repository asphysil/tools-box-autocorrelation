module time_correl
use nrtype 
use sys_constants, only: ndim 
use data_struc, only:  natms, cvel
USE nrtype; USE nrutil, ONLY :assert,assert_eq
use nr, only : cal_correl, correl_fft_dp

private
public :: cal_vel_correl, cal_phonon_dos, cal_pol_correl, cal_ir_spectra

contains

subroutine cal_vel_correl(nc, nvel, vv_correl)
implicit none 

integer, intent(in) :: nc , nvel  
real(dp), intent(out) :: vv_correl(nc) 

!local 
real(dp),dimension(:) :: cdata(nc), w(nc) 
real(dp) :: v0
integer :: i, j, k, l 
 

do i=1, nc 
     vv_correl(i) = 0.0 
     cdata(i) = 0.0 
 enddo 

do k =1, 3 
    do j=1, natms
      cdata(1:nvel) = cvel(1:nvel, j, k)
      cdata(nvel+1:nc) = 0.0
!            !print*, cdata(ncount)
!    enddo

     w = cal_correl(cdata, cdata)
     
       do l=1, nc 
        vv_correl(l) = vv_correl(l) + w(l)
       enddo
enddo 
enddo

!       !print*, 'AVG value ', vv_correl(1)
!       vv_correl(1:nc)= vv_correl(1:nc) !/vv_correl(1) ! normalization 
!     enddo
! enddo

v0 = vv_correl(1)
 do l=1, nc
     vv_correl(l) = vv_correl(l)/v0   ! normalization 
 enddo

    
end subroutine cal_vel_correl

subroutine cal_phonon_dos(nc, vv_corr, ph_dos)
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
 
!
  subroutine cal_pol_correl(nc, ndata, pp_corr)
    use data_struc, only : pol 
    implicit none 
    integer, intent(in) :: nc , ndata   
    real(dp), intent(out) :: pp_corr(nc) 
    
    !local 
    real(dp),dimension(:) :: cdata(nc), w(nc) 
    real(dp) :: v0
    integer :: i,  k
     
    
  
    pp_corr(1:) = 0.0 
    
    do k =1, 3 
         cdata(1:) =0.0
         cdata(1:ndata) = pol(k,1:ndata)
          
         w = cal_correl(cdata, cdata)
         
           do i=1, nc 
            pp_corr(i) = pp_corr(i) + w(i)
           enddo
    enddo 
    
    !       !print*, 'AVG value ', vv_correl(1)
    !       vv_correl(1:nc)= vv_correl(1:nc) !/vv_correl(1) ! normalization 
    !     enddo
    ! enddo
    
    v0 = pp_corr(1)
     do i=1, nc
         pp_corr(i) = pp_corr(i)/v0   ! normalization 
     enddo
    
        
    end subroutine cal_pol_correl  

    subroutine cal_ir_spectra(nc, pp_corr, ir_sp)
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
end module time_correl
