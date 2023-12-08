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
      cdata(1:nvel) = cvel(k, j, 1:nvel)
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

subroutine cal_phonon_dos(nvv_corr, vv_corr, nc, ph_dos)
    IMPLICIT NONE
    integer, intent(in) :: nc, nvv_corr  
    REAL(DP), DIMENSION(:), INTENT(IN) :: vv_corr(nvv_corr) 
    real(dp), dimension(:), intent(out) :: ph_dos(nc/2)

    COMPLEX(DPC), DIMENSION(nc/2) :: cdata 
    
real(dp), dimension(:) :: local_vv_corr(nc)
real(dp) :: n1_dp 
integer :: i 

local_vv_corr(1:nvv_corr) = vv_corr(1:nvv_corr)
local_vv_corr(nvv_corr+1:nc) = 0.0 

  cdata = correl_fft_dp(local_vv_corr) ! FFT 
 
  n1_dp = real(nvv_corr, kind=dp)

do i =1, nc/2 
    ph_dos(i) = abs(cdata(i))**2/n1_dp
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

    subroutine cal_ir_spectra(npp_corr, pp_corr, nc, ir_sp)
      IMPLICIT NONE
      integer, intent(in) :: nc, npp_corr 
      REAL(DP), DIMENSION(:), INTENT(IN) :: pp_corr(npp_corr) 
      real(dp), dimension(:), intent(out) :: ir_sp(nc/2)
  
      COMPLEX(DPC), DIMENSION(nc/2) :: cdata 
      
      real(dp), dimension(:) :: local_pp_corr(nc)
      integer :: i 
      real(dp) :: n1_dp 

      local_pp_corr(1:npp_corr) = pp_corr(1:npp_corr)
      local_pp_corr(npp_corr+1:nc) = 0.0 
  
    cdata = correl_fft_dp(local_pp_corr) ! FFT 
   
    n1_dp = real(npp_corr, kind=dp)

  do i =1, nc/2 
      ir_sp(i) = abs(cdata(i))**2/n1_dp
  enddo 
  
    end subroutine cal_ir_spectra   
end module time_correl
