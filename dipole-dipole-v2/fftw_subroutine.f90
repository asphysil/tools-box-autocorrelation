!=====================================
MODULE nrutil
  use nrtype 
  IMPLICIT NONE
INTEGER(I4B), PARAMETER :: NPAR_ARTH=16,NPAR2_ARTH=8 
!Each NPAR2 must be ≤ the corresponding NPAR.
INTEGER(I4B), PARAMETER :: NPAR_GEOP=4,NPAR2_GEOP=2
INTEGER(I4B), PARAMETER :: NPAR_CUMSUM=16
INTEGER(I4B), PARAMETER :: NPAR_CUMPROD=8
INTEGER(I4B), PARAMETER :: NPAR_POLY=8
INTEGER(I4B), PARAMETER :: NPAR_POLYTERM=8

INTERFACE swap
MODULE PROCEDURE swap_i,swap_r,swap_rv,swap_c, &
swap_cv,swap_cm,swap_z,swap_zv,swap_zm, &
masked_swap_rs,masked_swap_rv,masked_swap_rm
END INTERFACE

INTERFACE assert
MODULE PROCEDURE assert1,assert2,assert3,assert4,assert_v
END INTERFACE
INTERFACE assert_eq
MODULE PROCEDURE assert_eq2,assert_eq3,assert_eq4,assert_eqn
END INTERFACE

INTERFACE arth
MODULE PROCEDURE arth_r, arth_d, arth_i
END INTERFACE

contains 

SUBROUTINE swap_i(a,b)
  !Swap the contents of a and b.
  INTEGER(I4B), INTENT(INOUT) :: a,b
  INTEGER(I4B) :: dum
  dum=a
  a=b
  b=dum
END SUBROUTINE swap_i
SUBROUTINE swap_r(a,b)
    REAL(SP), INTENT(INOUT) :: a,b
    REAL(SP) :: dum
    dum=a
    a=b
    b=dum
END SUBROUTINE swap_r
SUBROUTINE swap_rv(a,b)
    REAL(SP), DIMENSION(:), INTENT(INOUT) :: a,b
    REAL(SP), DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
END SUBROUTINE swap_rv
SUBROUTINE swap_c(a,b)
    COMPLEX(SPC), INTENT(INOUT) :: a,b
    COMPLEX(SPC) :: dum
    dum=a
    a=b
    b=dum
END SUBROUTINE swap_c
SUBROUTINE swap_cv(a,b)
    COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: a,b
    COMPLEX(SPC), DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
END SUBROUTINE swap_cv
SUBROUTINE swap_cm(a,b)
    COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: a,b
    COMPLEX(SPC), DIMENSION(size(a,1),size(a,2)) :: dum
    dum=a
    a=b
    b=dum
END SUBROUTINE swap_cm
SUBROUTINE swap_z(a,b)
      COMPLEX(DPC), INTENT(INOUT) :: a,b
      COMPLEX(DPC) :: dum
      dum=a
      a=b
      b=dum
END SUBROUTINE swap_z
SUBROUTINE swap_zv(a,b)
      COMPLEX(DPC), DIMENSION(:), INTENT(INOUT) :: a,b
      COMPLEX(DPC), DIMENSION(SIZE(a)) :: dum
      dum=a
      a=b
      b=dum
END SUBROUTINE swap_zv
SUBROUTINE swap_zm(a,b)
      COMPLEX(DPC), DIMENSION(:,:), INTENT(INOUT) :: a,b
      COMPLEX(DPC), DIMENSION(size(a,1),size(a,2)) :: dum
      dum=a
      a=b
      b=dum
  END SUBROUTINE swap_zm
  SUBROUTINE masked_swap_rs(a,b,mask)
      REAL(SP), INTENT(INOUT) :: a,b
      LOGICAL(LGT), INTENT(IN) :: mask
      REAL(SP) :: swp
      if (mask) then
        swp=a
        a=b
        b=swp
        end if
  END SUBROUTINE masked_swap_rs
  SUBROUTINE masked_swap_rv(a,b,mask)
        REAL(SP), DIMENSION(:), INTENT(INOUT) :: a,b
        LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: mask
        REAL(SP), DIMENSION(size(a)) :: swp
        where (mask)
        swp=a
        a=b
        b=swp
        end where
  END SUBROUTINE masked_swap_rv
  SUBROUTINE masked_swap_rm(a,b,mask)
        REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a,b
        LOGICAL(LGT), DIMENSION(:,:), INTENT(IN) :: mask
        REAL(SP), DIMENSION(size(a,1),size(a,2)) :: swp
        where (mask)
        swp=a
        a=b
        b=swp
        end where
  END SUBROUTINE masked_swap_rm

! Routines for argument checking and error handling:
SUBROUTINE assert1(n1,string)
!Report and die if any logical is false (used for arg range checking).
CHARACTER(LEN=*), INTENT(IN) :: string
LOGICAL, INTENT(IN) :: n1
if (.not. n1) then
write (*,*) 'nrerror: an assertion failed with this tag:', &
string, '  program terminated by assert1'
STOP 
end if
END SUBROUTINE assert1

SUBROUTINE assert2(n1,n2,string)
  CHARACTER(LEN=*), INTENT(IN) :: string
  LOGICAL, INTENT(IN) :: n1,n2
  if (.not. (n1 .and. n2)) then
  write (*,*) 'nrerror: an assertion failed with this tag:', &
  string, 'program terminated by assert2'
  STOP 
  end if
  END SUBROUTINE assert2

  SUBROUTINE assert3(n1,n2,n3,string)
  CHARACTER(LEN=*), INTENT(IN) :: string
  LOGICAL, INTENT(IN) :: n1,n2,n3
  if (.not. (n1 .and. n2 .and. n3)) then
  write (*,*) 'nrerror: an assertion failed with this tag:', &
  string, 'program terminated by assert3'
  STOP 
  end if
  END SUBROUTINE assert3

  SUBROUTINE assert4(n1,n2,n3,n4,string)
  CHARACTER(LEN=*), INTENT(IN) :: string
  LOGICAL, INTENT(IN) :: n1,n2,n3,n4
  if (.not. (n1 .and. n2 .and. n3.and. n4)) then
  write (*,*) 'nrerror: an assertion failed with this tag:', &
  string,'program terminated by assert4'
  STOP 
  end if
  END SUBROUTINE assert4

SUBROUTINE assert_v(n,string)
CHARACTER(LEN=*), INTENT(IN) :: string
LOGICAL, DIMENSION(:), INTENT(IN) :: n
if (.not. all(n)) then
write (*,*) 'nrerror: an assertion failed with this tag:', &
string,'program terminated by assert_v'
STOP 
end if
END SUBROUTINE assert_v

FUNCTION assert_eq2(n1,n2,string)
  !Report and die if integers not all equal (used for size checking).
  CHARACTER(LEN=*), INTENT(IN) :: string
  INTEGER, INTENT(IN) :: n1,n2
  INTEGER :: assert_eq2
  if (n1 == n2) then
  assert_eq2=n1
  else
  write (*,*) 'nrerror: an assert_eq failed with this tag:', &
  string,'program terminated by assert_eq2'
  STOP 
  end if
END FUNCTION assert_eq2

FUNCTION assert_eq3(n1,n2,n3,string)
CHARACTER(LEN=*), INTENT(IN) :: string
INTEGER, INTENT(IN) :: n1,n2,n3
INTEGER :: assert_eq3
if (n1 == n2 .and. n2 == n3) then
assert_eq3=n1
else
write (*,*) 'nrerror: an assert_eq failed with this tag:', &
string,'program terminated by assert_eq3'
STOP 
end if
END FUNCTION assert_eq3

FUNCTION assert_eq4(n1,n2,n3,n4,string)
CHARACTER(LEN=*), INTENT(IN) :: string
INTEGER, INTENT(IN) :: n1,n2,n3,n4
INTEGER :: assert_eq4
if (n1 == n2 .and. n2 == n3.and. n3== n4) then
assert_eq4=n1
else
write (*,*) 'nrerror: an assert_eq failed with this tag:', &
string, 'program terminated by assert_eq4'
STOP 
end if
END FUNCTION assert_eq4

FUNCTION assert_eqn(nn,string)
CHARACTER(LEN=*), INTENT(IN) :: string
INTEGER, DIMENSION(:), INTENT(IN) :: nn
INTEGER :: assert_eqn
if (all(nn(2:) == nn(1))) then
assert_eqn=nn(1)
else
write (*,*) 'nrerror: an assert_eq failed with this tag:', &
string,'program terminated by assert_eqn'
STOP 
end if
END FUNCTION assert_eqn

FUNCTION zroots_unity(n,nn)
!Complex function returning nn powers of the nth root of unity.
  INTEGER(I4B), INTENT(IN) :: n,nn
  COMPLEX(SPC), DIMENSION(nn) :: zroots_unity
  INTEGER(I4B) :: k
  REAL(SP) :: theta
  zroots_unity(1)=1.0
  theta=TWOPI/n
  k=1
  do
  if (k >= nn) exit
  zroots_unity(k+1)=cmplx(cos(k*theta),sin(k*theta),SPC)
  zroots_unity(k+2:min(2*k,nn))=zroots_unity(k+1)*&
  zroots_unity(2:min(k,nn-k))
  k=2*k
  end do
  END FUNCTION zroots_unity

! Routines relating to polynomials and recurrences:
FUNCTION arth_r(first,increment,n)
!Array function returning an arithmetic progression.
REAL(SP), INTENT(IN) :: first,increment
INTEGER(I4B), INTENT(IN) :: n
REAL(SP), DIMENSION(n) :: arth_r
INTEGER(I4B) :: k,k2
REAL(SP) :: temp
if (n > 0) arth_r(1)=first
if (n <= NPAR_ARTH) then
do k=2,n
arth_r(k)=arth_r(k-1)+increment
end do
else
do k=2,NPAR2_ARTH
arth_r(k)=arth_r(k-1)+increment
end do
temp=increment*NPAR2_ARTH
k=NPAR2_ARTH
do
  if (k >= n) exit
  k2=k+k
  arth_r(k+1:min(k2,n))=temp+arth_r(1:min(k,n-k))
  temp=temp+temp
  k=k2
  end do
  end if
  END FUNCTION arth_r

  FUNCTION arth_d(first,increment,n)
    REAL(DP), INTENT(IN) :: first,increment
    INTEGER(I4B), INTENT(IN) :: n
    REAL(DP), DIMENSION(n) :: arth_d
    INTEGER(I4B) :: k,k2
    REAL(DP) :: temp
    if (n > 0) arth_d(1)=first
    if (n <= NPAR_ARTH) then
    do k=2,n
    arth_d(k)=arth_d(k-1)+increment
    end do
    else
    do k=2,NPAR2_ARTH
    arth_d(k)=arth_d(k-1)+increment
    end do
    temp=increment*NPAR2_ARTH
    k=NPAR2_ARTH
    do
    if (k >= n) exit
    k2=k+k
    arth_d(k+1:min(k2,n))=temp+arth_d(1:min(k,n-k))
    temp=temp+temp
    k=k2
    end do
    end if
    END FUNCTION arth_d

FUNCTION arth_i(first,increment,n)
INTEGER(I4B), INTENT(IN) :: first,increment,n
INTEGER(I4B), DIMENSION(n) :: arth_i
INTEGER(I4B) :: k,k2,temp
if (n > 0) arth_i(1)=first
if (n <= NPAR_ARTH) then
do k=2,n
arth_i(k)=arth_i(k-1)+increment
end do
else
do k=2,NPAR2_ARTH
arth_i(k)=arth_i(k-1)+increment
end do
temp=increment*NPAR2_ARTH
k=NPAR2_ARTH
do
if (k >= n) exit
k2=k+k
arth_i(k+1:min(k2,n))=temp+arth_i(1:min(k,n-k))
temp=temp+temp
k=k2
end do
end if
END FUNCTION arth_i

end module nrutil 

MODULE nr


  INTERFACE realft
  SUBROUTINE realft_dp(data,isign,zdata)
  USE nrtype
  REAL(DP), DIMENSION(:), INTENT(INOUT) :: data
  INTEGER(I4B), INTENT(IN) :: isign
  COMPLEX(DPC), DIMENSION(:), OPTIONAL, TARGET :: zdata
  END SUBROUTINE realft_dp

  SUBROUTINE realft_sp(data,isign,zdata)
  USE nrtype
  REAL(SP), DIMENSION(:), INTENT(INOUT) :: data
  INTEGER(I4B), INTENT(IN) :: isign
  COMPLEX(SPC), DIMENSION(:), OPTIONAL, TARGET :: zdata
  END SUBROUTINE realft_sp
  END INTERFACE

  INTERFACE fourrow
  SUBROUTINE fourrow_dp(data,isign)
    USE nrtype
    COMPLEX(DPC), DIMENSION(:,:), INTENT(INOUT) :: data
    INTEGER(I4B), INTENT(IN) :: isign
  END SUBROUTINE fourrow_dp

  SUBROUTINE fourrow_sp(data,isign)
    USE nrtype
    COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: data
    INTEGER(I4B), INTENT(IN) :: isign
  END SUBROUTINE fourrow_sp
END INTERFACE

  INTERFACE four1
    SUBROUTINE four1_dp(data,isign)
        USE nrtype
        COMPLEX(DPC), DIMENSION(:), INTENT(INOUT) :: data
        INTEGER(I4B), INTENT(IN) :: isign
    END SUBROUTINE four1_dp
    
    SUBROUTINE four1_sp(data,isign)
        USE nrtype
        COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: data
        INTEGER(I4B), INTENT(IN) :: isign
        END SUBROUTINE four1_sp
    END INTERFACE

  INTERFACE cal_correl
  FUNCTION correl_sp(data1,data2)
    USE nrtype
    REAL(SP), DIMENSION(:), INTENT(IN) :: data1,data2
    REAL(SP), DIMENSION(size(data1)) :: correl_sp 
    END FUNCTION correl_sp 

    FUNCTION correl_dp(data1,data2)
      USE nrtype
      REAL(DP), DIMENSION(:), INTENT(IN) :: data1,data2
      REAL(DP), DIMENSION(size(data1)) :: correl_dp
      END FUNCTION correl_dp 

  END INTERFACE cal_correl

  INTERFACE correl_fft_dp
    FUNCTION correl_fft_dp(data1)
      USE nrtype
      REAL(DP), DIMENSION(:), INTENT(IN) :: data1
      COMPLEX(DPC), DIMENSION(size(data1)/2) :: correl_fft_dp
      END FUNCTION correl_fft_dp 

  END INTERFACE correl_fft_dp

  INTERFACE twofft 
  SUBROUTINE twofft_sp(data1,data2, fft1, fft2)
    USE nrtype
    REAL(SP), DIMENSION(:), INTENT(IN) :: data1,data2
    COMPLEX(SP), DIMENSION(:), intent(out) :: fft1, fft2
    END SUBROUTINE twofft_sp 

    SUBROUTINE twofft_dp(data1,data2, fft1, fft2)
      USE nrtype
      REAL(DPC), DIMENSION(:), INTENT(IN) :: data1,data2
      COMPLEX(DPC), DIMENSION(:), intent(out) :: fft1, fft2
      END  SUBROUTINE twofft_dp

    END INTERFACE twofft
  end module nr 
  
!
!!! FFT
  SUBROUTINE fourrow_sp(data,isign)
    USE nrtype; USE nrutil, ONLY: assert,swap
    IMPLICIT NONE
    COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: data
    INTEGER(I4B), INTENT(IN) :: isign
    !Replaces each row (constant first index) of data(1:M ,1:N ) by its discrete Fourier trans-
    !form (transform on second index), if isign is input as 1; or replaces each row of data
    !by N times its inverse discrete Fourier transform, if isign is input as −1. N must be an
    ! integer power of 2. Parallelism is M -fold on the first index of data.
    INTEGER(I4B) :: n,i,istep,j,m,mmax,n2
    REAL(DP) :: theta
    COMPLEX(SPC), DIMENSION(size(data,1)) :: temp
    COMPLEX(DPC) :: w,wp !Double precision for the trigonometric recurrences.
    COMPLEX(SPC) :: ws
    n=size(data,2)
    call assert(iand(n,n-1)==0, 'n must be a power of 2 in fourrow_sp')
    n2=n/2
    j=n2
    !This is the bit-reversal section of the routine.
    do i=1,n-2
    if (j > i) call swap(data(:,j+1),data(:,i+1))
    m=n2
    do
    if (m < 2 .or. j < m) exit
    j=j-m
    m=m/2
    end do
    j=j+m
    end do
    mmax=1
    ! Here begins the Danielson-Lanczos section of the routine.
    ! do Outer loop executed log2 N times.
    do 
    if (n <= mmax) exit
    istep=2*mmax
    theta=PI_D/(isign*mmax) ! Initialize for the trigonometric recurrence.
    wp=cmplx(-2.0_dp*sin(0.5_dp*theta)**2,sin(theta),kind=dpc)
    w=cmplx(1.0_dp,0.0_dp,kind=dpc)

    do m=1,mmax !Here are the two nested inner loops.
      ws=w
      do i=m,n,istep
      j=i+mmax
      temp=ws*data(:,j) !This is the Danielson-Lanczos formula.
      data(:,j)=data(:,i)-temp
      data(:,i)=data(:,i)+temp
      end do
      w=w*wp+w !Trigonometric recurrence.
      end do
      mmax=istep
    end do
  END SUBROUTINE fourrow_sp

  SUBROUTINE fourrow_dp(data,isign)
    USE nrtype; USE nrutil, ONLY: assert,swap
    IMPLICIT NONE
    COMPLEX(DPC), DIMENSION(:,:), INTENT(INOUT) :: data
    INTEGER(I4B), INTENT(IN) :: isign
    INTEGER(I4B) :: n,i,istep,j,m,mmax,n2
    REAL(DP) :: theta
    COMPLEX(DPC), DIMENSION(size(data,1)) :: temp
    COMPLEX(DPC) :: w,wp
    COMPLEX(DPC) :: ws
    n=size(data,2)
    call assert(iand(n,n-1)==0, 'n must be a power of 2 in fourrow_dp')
    n2=n/2
    j=n2
    do i=1,n-2
    if (j > i) call swap(data(:,j+1),data(:,i+1))
    m=n2
    do
    if (m < 2 .or. j < m) exit
    j=j-m
    m=m/2
    end do
    j=j+m
    end do
    mmax=1

    do
    if (n <= mmax) exit
    istep=2*mmax
    theta=PI_D/(isign*mmax)
    wp=cmplx(-2.0_dp*sin(0.5_dp*theta)**2,sin(theta),kind=dpc)
    w=cmplx(1.0_dp,0.0_dp,kind=dpc)
    do m=1,mmax
    ws=w
    do i=m,n,istep
    j=i+mmax
    temp=ws*data(:,j)
    data(:,j)=data(:,i)-temp
    data(:,i)=data(:,i)+temp
    end do
    w=w*wp+w
    end do
    mmax=istep
    end do
  END SUBROUTINE fourrow_dp

  SUBROUTINE four1_sp(data,isign)
    USE nrtype
    USE nrutil, ONLY: arth,assert
    USE nr, ONLY: fourrow
    IMPLICIT NONE
    COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(I4B), INTENT(IN) :: isign
    ! Replaces a complex array data by its discrete Fourier transform, if isign is input as 1;
    ! or replaces data by its inverse discrete Fourier transform times the size of data, if isign
    ! is input as −1. The size of data must be an integer power of 2. Parallelism is achieved
    ! by internally reshaping the input array to two dimensions. (Use this version if fourrow is
    ! faster than fourcol on your machine.)
    COMPLEX(SPC), DIMENSION(:,:), ALLOCATABLE :: dat,temp
    COMPLEX(DPC), DIMENSION(:), ALLOCATABLE :: w,wp
    REAL(DP), DIMENSION(:), ALLOCATABLE :: theta
    INTEGER(I4B) :: n,m1,m2,j
    n=size(data)
    call assert(iand(n,n-1)==0, 'n must be a power of 2 in four1_sp')
    !Find dimensions as close to square as possible, allocate space, and reshape the input array.
    m1=2**ceiling(0.5_sp*log(real(n,sp))/0.693147_sp)
    m2=n/m1
    allocate(dat(m1,m2),theta(m1),w(m1),wp(m1),temp(m2,m1))
    dat=reshape(data,shape(dat))
    call fourrow(dat,isign) 
    theta=arth(0,isign,m1)*TWOPI_D/n !Set up recurrence.
    wp=cmplx(-2.0_dp*sin(0.5_dp*theta)**2,sin(theta),kind=dpc)
    w=cmplx(1.0_dp,0.0_dp,kind=dpc)
    do j=2,m2 !Multiply by the extra phase factor.
    w=w*wp+w
    dat(:,j)=dat(:,j)*w
    end do
    temp=transpose(dat) !Transpose, and transform on (original) first index.
    call fourrow(temp,isign)
    data=reshape(temp,shape(data))! Reshape the result back to one dimension.
    deallocate(dat,w,wp,theta,temp)
    END SUBROUTINE four1_sp

    SUBROUTINE four1_dp(data,isign)
      USE nrtype; USE nrutil, ONLY: arth,assert
      USE nr, ONLY: fourrow
      IMPLICIT NONE
      COMPLEX(DPC), DIMENSION(:), INTENT(INOUT) :: data
      INTEGER(I4B), INTENT(IN) :: isign
      COMPLEX(DPC), DIMENSION(:,:), ALLOCATABLE :: dat,temp
      COMPLEX(DPC), DIMENSION(:), ALLOCATABLE :: w,wp
      REAL(DP), DIMENSION(:), ALLOCATABLE :: theta
      INTEGER(I4B) :: n,m1,m2,j
      n=size(data)
      call assert(iand(n,n-1)==0, 'n must be a power of 2 in four1_dp')
      m1=2**ceiling(0.5_sp*log(real(n,sp))/0.693147_sp)
      m2=n/m1
      allocate(dat(m1,m2),theta(m1),w(m1),wp(m1),temp(m2,m1))
      dat=reshape(data,shape(dat))
      call fourrow(dat,isign)
      theta=arth(0,isign,m1)*TWOPI_D/n
      wp=cmplx(-2.0_dp*sin(0.5_dp*theta)**2,sin(theta),kind=dpc)
      w=cmplx(1.0_dp,0.0_dp,kind=dpc)
      do j=2,m2
      w=w*wp+w
      dat(:,j)=dat(:,j)*w
      end do
      temp=transpose(dat)
      call fourrow(temp,isign)
      data=reshape(temp,shape(data))
      deallocate(dat,w,wp,theta,temp)
      END SUBROUTINE four1_dp



  SUBROUTINE realft_sp(data,isign,zdata)
    USE nrtype
    USE nrutil, ONLY: assert,assert_eq,zroots_unity
    USE nr, ONLY: four1
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(I4B), INTENT(IN) :: isign
    COMPLEX(SPC), DIMENSION(:), OPTIONAL, TARGET :: zdata
    ! When isign = 1, calculates the Fourier transform of a set of N real-valued data points,
    ! input in the array data. If the optional argument zdata is not present, the data are replaced
    ! by the positive frequency half of its complex Fourier transform. The real-valued first and
    ! last components of the complex transform are returned as elements data(1) and data(2),
    ! respectively. If the complex array zdata of length N/2 is present, data is unchanged and
    ! the transform is returned in zdata. N must be a power of 2. If isign = −1, this routine
    ! calculates the inverse transform of a complex data array if it is the transform of real data.
    ! (Result in this case must be multiplied by 2/N .) The data can be supplied either in data,
    ! with zdata absent, or in zdata.
    INTEGER(I4B) :: n,ndum,nh,nq
    COMPLEX(SPC), DIMENSION(size(data)/4) :: w
    COMPLEX(SPC), DIMENSION(size(data)/4-1) :: h1,h2
    COMPLEX(SPC), DIMENSION(:), POINTER :: cdata 
    ! Used for internal complex computations.
    
    COMPLEX(SPC) :: z
    REAL(SP) :: c1=0.5_sp,c2
    n=size(data)
    call assert(iand(n,n-1)==0, 'n must be a power of 2 in realft_sp')
    nh=n/2
    nq=n/4
    if (present(zdata)) then
      ndum=assert_eq(n/2,size(zdata),'realft_sp')
      cdata=>zdata !Use zdata as cdata.
      if (isign == 1) cdata=cmplx(data(1:n-1:2),data(2:n:2),kind=spc)
      else
      allocate(cdata(n/2)) !Have to allocate storage ourselves.
      cdata=cmplx(data(1:n-1:2),data(2:n:2),kind=spc)
      end if
      if (isign == 1) then
      c2=-0.5_sp
      call four1(cdata,+1) 
      !The forward transform is here.
      ! Otherwise set up for an inverse trans-form.
      else 
        
        c2=0.5_sp
      end if
      w=zroots_unity(sign(n,isign),n/4)
      w=cmplx(-aimag(w),real(w),kind=spc)
      h1=c1*(cdata(2:nq)+conjg(cdata(nh:nq+2:-1))) 
      !The two separate transforms are sep-arated out of cdata.
      h2=c2*(cdata(2:nq)-conjg(cdata(nh:nq+2:-1)))
      !Next they are recombined to form the true transform of the original real data:
      cdata(2:nq)=h1+w(2:nq)*h2
      cdata(nh:nq+2:-1)=conjg(h1-w(2:nq)*h2)
      z=cdata(1) 
      !Squeeze the first and last data to-gether to get them all within the
      ! original array.
      if (isign == 1) then
      cdata(1)=cmplx(real(z)+aimag(z),real(z)-aimag(z),kind=spc)
      else
      cdata(1)=cmplx(c1*(real(z)+aimag(z)),c1*(real(z)-aimag(z)),kind=spc)
      call four1(cdata,-1) 
      !This is the inverse transform for the case isign=-1
      end if
      if (present(zdata)) then 
        ! Ship out answer in data if required.
      if (isign /= 1) then
      data(1:n-1:2)=real(cdata)
      data(2:n:2)=aimag(cdata)
      end if
      else
      data(1:n-1:2)=real(cdata)
      data(2:n:2)=aimag(cdata)
      deallocate(cdata)
      end if
    END SUBROUTINE realft_sp

  SUBROUTINE realft_dp(data,isign,zdata)
    USE nrtype; USE nrutil, ONLY: assert,assert_eq,zroots_unity
    USE nr, ONLY: four1
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(I4B), INTENT(IN) :: isign
    COMPLEX(DPC), DIMENSION(:), OPTIONAL, TARGET :: zdata
    INTEGER(I4B) :: n,ndum,nh,nq
    COMPLEX(DPC), DIMENSION(size(data)/4) :: w
    COMPLEX(DPC), DIMENSION(size(data)/4-1) :: h1,h2
    COMPLEX(DPC), DIMENSION(:), POINTER :: cdata
    COMPLEX(DPC) :: z
    REAL(DP) :: c1=0.5_dp,c2
    n=size(data)
    call assert(iand(n,n-1)==0, 'n must be a power of 2 in realft_dp')
    nh=n/2
    nq=n/4
    if (present(zdata)) then
    ndum=assert_eq(n/2,size(zdata),'realft_dp')
    cdata=>zdata
    if (isign == 1) cdata=cmplx(data(1:n-1:2),data(2:n:2),kind=spc)
    else
    allocate(cdata(n/2))
    cdata=cmplx(data(1:n-1:2),data(2:n:2),kind=spc)
    end if
    if (isign == 1) then
    c2=-0.5_dp
    call four1(cdata,+1)
    else
    c2=0.5_dp
    end if
    w=zroots_unity(sign(n,isign),n/4)
    w=cmplx(-aimag(w),real(w),kind=dpc)
    h1=c1*(cdata(2:nq)+conjg(cdata(nh:nq+2:-1)))
    h2=c2*(cdata(2:nq)-conjg(cdata(nh:nq+2:-1)))
    cdata(2:nq)=h1+w(2:nq)*h2
    cdata(nh:nq+2:-1)=conjg(h1-w(2:nq)*h2)
    z=cdata(1)
    if (isign == 1) then
    cdata(1)=cmplx(real(z)+aimag(z),real(z)-aimag(z),kind=dpc)
    else
      cdata(1)=cmplx(c1*(real(z)+aimag(z)),c1*(real(z)-aimag(z)),kind=dpc)
      call four1(cdata,-1)
      end if
      if (present(zdata)) then
      if (isign /= 1) then
      data(1:n-1:2)=real(cdata)
      data(2:n:2)=aimag(cdata)
      end if
      else
      data(1:n-1:2)=real(cdata)
      data(2:n:2)=aimag(cdata)
      deallocate(cdata)
      end if
  END SUBROUTINE realft_dp

  FUNCTION correl_sp(data1,data2)
    USE nrtype; USE nrutil, ONLY :assert,assert_eq
    USE nr, ONLY :realft
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(IN) :: data1,data2

    real(SP), dimension(size(data1)) :: loc_data1, loc_data2

    REAL(SP), DIMENSION(size(data1)) :: correl_sp
   ! Computes the correlation of two real data sets data1 and data2 of length N (includ-
   ! ing any user-supplied zero padding). N must be an integer power of 2. The answer is
   ! returned as the function correl, an array of length N . The answer is stored in wrap-
   ! around order, i.e., correlations at increasingly negative lags are in correl(N ) on down to
   ! correl(N/2 + 1), while correlations at increasingly positive lags are in correl(1) (zero
   ! lag) on up to correl(N/2). Sign convention of this routine: if data1 lags data2, i.e.,
   ! is shifted to the right of it, then correl will show a peak at positive lags.
    COMPLEX(SPC), DIMENSION(size(data1)/2) :: cdat1,cdat2
    INTEGER(I4B) :: no2,n !Normalization for inverse FFT.
    n=assert_eq(size(data1),size(data2),'correl_sp')
    call assert(iand(n,n-1)==0, 'n must be a power of 2 in correl_sp')
    no2=n/2

    loc_data1(1:) = data1(1:)
    loc_data2(1:) = data2(1:)

    call realft(loc_data1,1,cdat1) !Transform both data vectors.
    call realft(loc_data2,1,cdat2)
    cdat1(1)=cmplx(real(cdat1(1))*real(cdat2(1))/no2, & 
             aimag(cdat1(1))*aimag(cdat2(1))/no2, kind=spc)
    !  Multiply to find FFT of their
    cdat1(2:)=cdat1(2:)*conjg(cdat2(2:))/no2
    call realft(correl_sp,-1,cdat1) ! Inverse transform gives correlation.
  END FUNCTION correl_sp 

  FUNCTION correl_dp(data1,data2)
    USE nrtype; USE nrutil, ONLY :assert,assert_eq
    USE nr, ONLY :realft
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN) :: data1,data2
    real(DP), dimension(size(data1)) :: loc_data1, loc_data2
    REAL(DP), DIMENSION(size(data1)) :: correl_dp
   ! Computes the correlation of two real data sets data1 and data2 of length N (includ-
   ! ing any user-supplied zero padding). N must be an integer power of 2. The answer is
   ! returned as the function correl, an array of length N . The answer is stored in wrap-
   ! around order, i.e., correlations at increasingly negative lags are in correl(N ) on down to
   ! correl(N/2 + 1), while correlations at increasingly positive lags are in correl(1) (zero
   ! lag) on up to correl(N/2). Sign convention of this routine: if data1 lags data2, i.e.,
   ! is shifted to the right of it, then correl will show a peak at positive lags.
    COMPLEX(DPC), DIMENSION(size(data1)/2) :: cdat1,cdat2
    INTEGER(I4B) :: no2,n !Normalization for inverse FFT.
    n=assert_eq(size(data1),size(data2),'correl_dp')
    call assert(iand(n,n-1)==0, 'n must be a power of 2 in correl_dp')
    no2=n/2
        
    loc_data1(1:) = data1(1:)
    loc_data2(1:) = data2(1:)

    call realft(loc_data1,1,cdat1) !Transform both data vectors
    call realft(loc_data2,1,cdat2)

    cdat1(1)=cmplx(real(cdat1(1))*real(cdat2(1))/no2, & 
             aimag(cdat1(1))*aimag(cdat2(1))/no2, kind=DPC)
    !  Multiply to find FFT of their
    cdat1(2:)=cdat1(2:)*conjg(cdat2(2:))/no2

    call realft(correl_dp,-1,cdat1) ! Inverse transform gives correlation.
  END FUNCTION correl_dp 




  FUNCTION correl_fft_dp(data1)
    USE nrtype; USE nrutil, ONLY :assert,assert_eq
    USE nr, ONLY :realft
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN) :: data1
    real(DP), dimension(size(data1)) :: loc_data1
    COMPLEX(DPC), DIMENSION(size(data1)/2) :: correl_fft_dp
   ! Computes the correlation of two real data sets data1 and data2 of length N (includ-
   ! ing any user-supplied zero padding). N must be an integer power of 2. The answer is
   ! returned as the function correl, an array of length N . The answer is stored in wrap-
   ! around order, i.e., correlations at increasingly negative lags are in correl(N ) on down to
   ! correl(N/2 + 1), while correlations at increasingly positive lags are in correl(1) (zero
   ! lag) on up to correl(N/2). Sign convention of this routine: if data1 lags data2, i.e.,
   ! is shifted to the right of it, then correl will show a peak at positive lags.

    INTEGER(I4B) :: no2,n !Normalization for inverse FFT.
    n=assert_eq(size(data1),size(data1),'correl_dp')
    call assert(iand(n,n-1)==0, 'n must be a power of 2 in correl_dp')
    no2=n/2
  
    loc_data1(1:) = data1(1:)


    call realft(loc_data1,1,correl_fft_dp) !Transform both data vectors


   ! correl_fft_dp(1)=cmplx(real(correl_fft_dp(1))/no2, & 
   !         aimag(correl_fft_dp(1))/no2, kind=DPC)
    !  Multiply to find FFT of their
   ! correl_fft_dp(2:)=correl_fft_dp(2:)/no2

    correl_fft_dp(1)=cmplx(real(correl_fft_dp(1)), & 
            aimag(correl_fft_dp(1)), kind=DPC)
    !  Multiply to find FFT of their
    correl_fft_dp(2:)=correl_fft_dp(2:)
  END FUNCTION correl_fft_dp 


  SUBROUTINE twofft_sp(data1,data2,fft1,fft2)
    USE nrtype; USE nrutil, ONLY: assert,assert_eq
    USE nr, ONLY: four1
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(IN) :: data1,data2
    COMPLEX(SPC), DIMENSION(:), INTENT(OUT) :: fft1,fft2
    ! Given two real input arrays data1 and data2 of length N , this routine calls four1 and
    ! returns two complex output arrays, fft1 and fft2, each of complex length N , that contain
    ! the discrete Fourier transforms of the respective data arrays. N must be an integer power
    ! of 2.
    INTEGER(I4B) :: n,n2
    COMPLEX(SPC), PARAMETER :: C1=(0.5_sp,0.0_sp), C2=(0.0_sp,-0.5_sp)
    COMPLEX(SPC), DIMENSION(size(data1)/2+1) :: h1,h2

    n=assert_eq(size(data1),size(data2),size(fft1),size(fft2),'twofft_sp')
    call assert(iand(n,n-1)==0, 'n must be a power of 2 in twofft_sp')
    fft1=cmplx(data1,data2,kind=spc) !Pack the two real arrays into one complex array.
    call four1(fft1,1) !Transform the complex array.
    fft2(1)=cmplx(aimag(fft1(1)),0.0_sp,kind=spc)
    fft1(1)=cmplx(real(fft1(1)),0.0_sp,kind=spc)
    n2=n/2+1
    h1(2:n2)=C1*(fft1(2:n2)+conjg(fft1(n:n2:-1))) !Use symmetries to separate the two transforms.
    h2(2:n2)=C2*(fft1(2:n2)-conjg(fft1(n:n2:-1)))
    fft1(2:n2)=h1(2:n2) ! Ship them out in two complex arrays.
    fft1(n:n2:-1)=conjg(h1(2:n2))
    fft2(2:n2)=h2(2:n2)
    fft2(n:n2:-1)=conjg(h2(2:n2))
    END SUBROUTINE twofft_sp

    SUBROUTINE twofft_dp(data1,data2,fft1,fft2)
      USE nrtype; USE nrutil, ONLY: assert,assert_eq
      USE nr, ONLY: four1
      IMPLICIT NONE
      REAL(DP), DIMENSION(:), INTENT(IN) :: data1,data2
      COMPLEX(DPC), DIMENSION(:), INTENT(OUT) :: fft1,fft2
      ! Given two real input arrays data1 and data2 of length N , this routine calls four1 and
      ! returns two complex output arrays, fft1 and fft2, each of complex length N , that contain
      ! the discrete Fourier transforms of the respective data arrays. N must be an integer power
      ! of 2.
      INTEGER(I4B) :: n,n2
      COMPLEX(DPC), PARAMETER :: C1=(0.5_sp,0.0_sp), C2=(0.0_sp,-0.5_sp)
      COMPLEX(DPC), DIMENSION(size(data1)/2+1) :: h1, h2
      n=assert_eq(size(data1),size(data2),size(fft1),size(fft2),'twofft_sp')
      call assert(iand(n,n-1)==0, 'n must be a power of 2 in twofft_sp')
      fft1=cmplx(data1,data2,kind=spc) !Pack the two real arrays into one complex array.
      call four1(fft1,1) !Transform the complex array.
      fft2(1)=cmplx(aimag(fft1(1)),0.0_sp,kind=DPC)
      fft1(1)=cmplx(real(fft1(1)),0.0_sp,kind=DPC)
      n2=n/2+1
      h1(2:n2)=C1*(fft1(2:n2)+conjg(fft1(n:n2:-1))) !Use symmetries to separate the two transforms.
      h2(2:n2)=C2*(fft1(2:n2)-conjg(fft1(n:n2:-1)))
      fft1(2:n2)=h1(2:n2) ! Ship them out in two complex arrays.
      fft1(n:n2:-1)=conjg(h1(2:n2))
      fft2(2:n2)=h2(2:n2)
      fft2(n:n2:-1)=conjg(h2(2:n2))
      END SUBROUTINE twofft_dp 
