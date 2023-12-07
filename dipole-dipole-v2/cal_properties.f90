module sys_properties_cal 
    use nrtype 
    use sys_constants, only: ndim, Navogadro, ev2J, kb   
    use data_struc, only: dt, atms_mass, latt_vec, pos
    implicit none
    private 
    public :: cal_velocity_forward_method, cal_polarization
contains

subroutine cal_velocity_forward_method(natms, nskip, nvel)
    use data_struc, only :  atms_vel, cvel, vel_fdiff, ek, temp 
    implicit none 
    integer, intent(in) :: natms, nskip, nvel  

! local 
integer :: i, j, k
real(dp) :: idt 
real(dp) :: x, n_dp  
real(dp) :: v(3)
real(dp) :: factor1

! Calculating velocity using forward difference 
idt = 1.0/dt 
! can be made more efficient this section
do i=1, nvel 
    do j =1, natms 
       do k =1, 3 
        v(k)= (pos(k, j, nskip+i+1)- pos(k, j, nskip+i))
        ! PBC applied to account correct velocity direction 
        if (v(k)>0.5 ) then  ! 
            v(k)= v(k) - 1.0 
        elseif (v(k) <-0.5 ) then 
            v(k) = v(k) + 1.0 
        endif 

       enddo 
       vel_fdiff(1:3,j, i) = matmul(v, latt_vec) *idt ! fraction to xyz
       atms_vel(1:3,j,i) = vel_fdiff(1:3,j,i)
    enddo
 
enddo

! calculating K.E and T
factor1 = real(natms, kind=dp) ! N 
do i= 1, nvel
    x=0.0 
    do j=1, natms
        v(1:3) = atms_vel(1:3, j, i)
        x = x + 0.5*atms_mass(j)*dot_product(v,v) ! m*v^2/2 
    enddo 
!  (Angs/fs)**2 = (m/s)**2 * 1.0E10
!  1amu = 1/Navogadro = gm = (kg)*1.0E-3
! amu* (Angs/fs)**2 = (kg)*1.0E-3*(m/s)**2 * 1.0E10 = 1.0E7 Joule 
    !
    ek(i) = (x *1.0E7)/Navogadro ! Joule
    ek(i) = ek(i)/ev2J ! in eV 
    temp(i) = ek(i)*2.0/(3.0*kb*factor1)! in K
enddo 

! substracting average velocity 
n_dp = real(nvel, kind=dp)

do k= 1, 3
    do j=1, natms
        x = sum(vel_fdiff(k,j,1:), dim=1)/n_dp ! Average velocity
        !print*,x 
        do i = 1, nvel 
            cvel(k,j,i) = vel_fdiff(k,j,i) - x 
        enddo
    enddo
enddo

end subroutine cal_velocity_forward_method

 subroutine cal_polarization(natms, nskip, ndata)
    use data_struc, only : pol, born_charg 

    integer, intent(in) :: natms, nskip, ndata 


!local
    integer :: i, j, k 
real(dp) :: v_frac(3), v_car(3)
real(dp) :: dpol(3), tot_dpol(3)
real(dp) :: z_star(3,3)

real(dp) :: np_dp 
real(dp) :: avg_p(3) 

real(dp) :: na_dp1 

    do i=1, ndata
        tot_dpol(1:3) = 0.0 

        ! can be make eficient way to calculate polarization 
        do j =1, natms  
            do k =1, 3
                v_frac(k) = pos(k, j, nskip+i) 
                 ! PBC applied to account correct velocity direction 
                if (v_frac(k)>0.5 ) then  ! 
                   v_frac(k)= v_frac(k) - 1.0 
                elseif (v_frac(k) <-0.5 ) then 
                    v_frac(k) = v_frac(k) + 1.0 
                endif 
            enddo
                v_car(1:3) = matmul(v_frac, latt_vec)
                z_star = born_charg(:,:,j) 
                dpol = matmul(v_car, z_star)
                tot_dpol(1:3) = tot_dpol(1:3) + dpol(1:3)
        enddo
        pol(1:3, i) = tot_dpol(1:3)
    enddo

    np_dp = real(ndata, kind=dp)

    do k=1, 3
        avg_p(k) = sum(pol(k,:))/ np_dp 
    enddo

        do i = 1, ndata 
            do k =1 ,3 
                  pol(k,i) = pol(k,i) - avg_p(k)
            enddo
        enddo

 end subroutine cal_polarization   

end module sys_properties_cal 
