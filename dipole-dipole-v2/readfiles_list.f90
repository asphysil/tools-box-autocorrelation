MODULE readfiles
 use nrtype
 use sys_constants, only : ndim 
 use data_struc, only : ntype_atms, latt_vec
    IMPLICIT NONE 
private
public :: read_vasp_xdatcar,read_test_data

contains

subroutine read_vasp_xdatcar(fileid, ntype,  natms, ndata,  atms_pos)
   implicit none 
   integer, intent(in) :: fileid, ndata, ntype, natms 
   real(dp), intent(out) :: atms_pos(ndim, natms, ndata)   

   ! local 
   real(dp) :: sl 
   real(dp) :: v(3)
   integer :: loc_natms 
   integer :: i, j, k
   character(len=2) :: atm_name 

print*, 'Reading file ...'

   read(fileid, *)
   read(fileid, *) sl 

   do i = 1, 3
      read(fileid, *) (latt_vec(i, j), j=1,3)
   enddo

   do i = 1, 3
      do j = 1, 3 
         latt_vec(i, j) = latt_vec(i,j)*sl 
      enddo
   enddo 
   
print*, ' Lattice parametrs '
do i=1,3
print*, (latt_vec(i, j), j=1,3)
enddo

   read(fileid, *) (atm_name, i=1, ntype)
   !print*, " element name  ", i, ' = ', atm_name, ' out of ', ntype
   !enddo 

   read(fileid, *) (ntype_atms(j), j=1, ntype)
   !print*, ( ' Element No.', j, ' number of atom ' , ntype_atms(j), ' out of ', ntype,  j=1, ntype)
    print*, (ntype_atms(j),j=1, ntype)

   loc_natms = sum(ntype_atms, dim=1)

if (loc_natms .ne. natms) then 
    print*, ' wrong total number atoms in the inputs and ouptput  '
else
   print*, "total number atoms in the inputs and ouptput correct"
endif 

do i = 1, ndata
    
  read(fileid, *)
   do j=1, natms 
     read(fileid, *) (atms_pos(k, j, i), k=1, 3)
     !(v(k), k=1,3)
    ! atms_pos(1:3, j, i) = v(1:3)
     !matmul(v, latt_vec)
     !print*, (atms_pos(k, j, i), k=1, 3)
   enddo 
enddo 
print*, '----End reading ---'
end subroutine read_vasp_xdatcar

subroutine read_test_data(fileID,  ndata, x, y)
   integer, intent(in) :: fileID, ndata
   real(dp), intent(out) ::  x(ndata), y(ndata)

   integer :: i 

   do i=1, ndata
      read(fileID, *) x(i), y(i) 
   enddo
end subroutine read_test_data
END MODULE readfiles
