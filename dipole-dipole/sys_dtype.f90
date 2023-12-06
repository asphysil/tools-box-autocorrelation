module sys_constants
    use nrtype 
    implicit none
    integer, parameter :: ndim = 3  
    !real(dp), parameter :: t_fs2sec=1.0e-15 ! fs to sec
    !real(dp), parameter :: angs2m = 1.0e-10 ! Angstrom to meter 
    !real(dp), parameter :: 1.660539e-27 ! amu to kg 
    ! v = x/t = 1.0e-10/1.0e-15 = 1.0e5
    ! m*v**2 =1.660539e-27*1.0e10=1.660539e-17  ! jule
    ! jule2eV=6.241509074461e+18
    ! 1.660539e-17*6.241509074461e+18=6.241509074461*16.60539
    real(dp), parameter :: fs2sec=1.0E-15 
    ! # Boltzmann Constant in [eV/K]
    real(dp), parameter :: kb = 8.617332478E-5
    ! electron volt in [Joule]
    real(dp), parameter :: ev2J = 1.60217733E-19 ! eV2J
    ! Avogadro's Constant
    real(dp), parameter :: Navogadro = 6.0221412927E23
end module sys_constants

module filesystem 
    implicit none
    integer :: read_fileid 
    integer :: write_fileid1
    integer :: write_fileid2
    integer :: write_fileid3
    integer :: write_fileid4

end module filesystem 

module data_struc
    use nrtype 
    use sys_constants, only :  ndim 
    implicit none

    !

    !
    integer :: ntype, natms  
    integer :: ntot, ndata, nvdata, nc, nc_fft   
    
    integer :: nc_ir, nc_ir_fft ! for ir_spectra 

    integer, dimension(:), allocatable :: ntype_atms

    real(dp) :: dt 

    real(dp), dimension(:,:), allocatable :: latt_vec 
    real(dp), dimension(:,:,:), allocatable :: pos, atms_vel 
    real(dp), dimension(:,:,:), allocatable ::  vel_cdiff, vel_fdiff !,vel_bdiff
    real(dp), dimension(:, :, :), allocatable :: cvel 
  
    real(dp), dimension(:), allocatable :: vv_corr

    real(DP),dimension(:), allocatable :: ph_dos

    ! For fft test 
    real(dp), dimension(:), allocatable :: atms_mass 
    real(dp), dimension(:), allocatable :: ek, temp 
    real(dp), dimension(:,:), allocatable :: pol 
    real(dp), dimension(:, :, :), allocatable :: born_charg
    real(dp), dimension(:), allocatable :: pp_corr
    real(dp), dimension(:), allocatable :: ir_sp 
end module data_struc

! module test_data 
!     use sys_constants, only: ndim 
!     use nrtype
!     implicit none
!     integer :: ntest 
    
!     real(dp) :: dt_test 
!     real(dp), dimension(:), allocatable :: t_test,  y_test, x_test, w_test 
! end module test_data 