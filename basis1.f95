! basis creation
!
module basis
use types
use constants, only: I, pi
use strings
use fact
use dispmodule
implicit none

type Hbas ! Hermitian basis
    complex(dp), dimension(:,:), allocatable :: matrix
    character(8) idx
    logical :: hermitian = .true.
    logical :: unitary = .false.
end type Hbas

type Ubas ! Unitary basis
    complex(dp), dimension(:,:), allocatable :: matrix
    character(8) idx
    logical :: unitary = .true.
    logical :: hermitian = .false.
end type Ubas

type HBasis
    type(Hbas), dimension(:), allocatable :: bas
    integer dimen
end type HBasis

type UBasis
    type(Ubas), dimension(:), allocatable :: bas
    integer dimen
end type UBasis

type comp
    character(6) idx
    real :: x
end type comp

type components
    type(comp), dimension(3) :: cmp
end type components

interface printBasis
    module procedure printHbasis
    module procedure printUbasis
end interface

contains

function E(j,k,d) result(rr)
integer, intent(in) :: j,k,d
complex(dp), dimension(:,:), allocatable :: rr

allocate(rr(d,d))

rr=0
rr(j,k)=1

end function E

function f(k,j,d) result(rr)
integer, intent(in) :: j,k,d
complex(dp), dimension(:,:), allocatable :: rr

allocate(rr(d,d))

rr=0
if (k < j) then
    rr=E(k,j,d)+E(j,k,d)
else if (k > j) then
    rr=-I*(E(j,k,d)-E(k,j,d))
endif

end function f

function id(d) result(rr)
integer,intent(in) :: d
integer :: i
complex(dp), dimension(:,:), allocatable :: rr

allocate(rr(d,d))

rr=0
forall (i=1:d) rr(i,i)=1

end function id

function dirsum(m1,s2) result(rr)
complex(dp), dimension(:,:), intent(in) :: m1
integer, intent(in) :: s2
complex(dp), dimension(:,:), allocatable :: rr
integer, dimension(2) :: shap, sh1
integer :: i,j

sh1=shape(m1)
shap=sh1+1
allocate(rr(shap(1),shap(2)))

rr=0
forall (i=1:sh1(1), j=1:sh1(2)) rr(i,j)=m1(i,j)
rr(sh1(1)+1,sh1(2)+1)=s2

end function dirsum

recursive function h(k,d) result(rr)
integer, intent(in) :: k,d
complex(dp), dimension(:,:), allocatable :: rr
real :: c1

c1=sqrt(2.0/(d*(d-1)))

if (k == 1) then
    rr=id(d)
else if (k == d) then
    rr=c1*(dirsum(h(1,d-1),1-d))
else
    rr=dirsum(h(k,d-1),0)
endif

end function h

function clockm(d) result(rr)
    integer, intent(in) :: d
    complex(dp), dimension(d,d) :: rr
    complex(dp) :: omega,ompar
    integer :: i

    ompar=complex(0d0,2.0d0*pi/d)
    omega=zexp(ompar)

    rr=0
    forall (i=0:d-1) rr(i+1,i+1)=omega**i
end function clockm

function skewm(d) result(rr)
    integer, intent(in) :: d
    complex(dp), dimension(d,d) :: rr
    integer :: i

    rr=0
    rr(1,d)=1
    forall (i=1:d-1) rr(i+1,i)=1
end function skewm

function genPauli(j,k,d) result(rr)
    integer, intent(in) :: j,k,d
    integer :: i
    complex(dp), dimension(d,d) :: rr,cm,sm,t1,t2,t3,t4

    cm=clockm(d)
    sm=skewm(d)

    t1=cm
    t2=sm
    t3=cm
    t4=sm

    if ((j==1) .and. (k==1)) then
        forall(i=1:d) rr(i,i)=1
    else if (j==1) then
        do i=1,d
            t3=matmul(t3,cm)
        end do
        rr=t3
    else if (k==1) then
        do i=1,d
            t4=matmul(t4,sm)
        end do
        rr=t4
    else if ((k>1) .and. (j>1)) then
        do i=1,k
            t1=matmul(t1,cm)
        end do
        do i=1,j
            t2=matmul(t2,sm)
        end do
        rr=matmul(t1,t2)
    end if
end function genPauli

function genHbasis(d) result(rrbas)
    integer, intent(in) :: d
    type(Hbas) :: b1
    type(Hbasis) :: rrbas
    integer :: i,j,counti,ed
    character(len=2) c,c1,c2

    ed=perm(d,2)+d
    allocate(rrbas%bas(ed))

    print *,'Generating basis ...'

    do i=1,ed
        allocate(rrbas%bas(i)%matrix(d,d))
    end do

    call writenum(d,c1,'i2')
    call removesp(c1)

! H matrices
    do i=1,d
        b1%matrix=h(i,d)
        call writenum(i,c,'i2')
        call removesp(c)
        b1%idx='h('//c//','//c1//')'
        call removesp(b1%idx)
        rrbas%bas(i)=b1
    end do
! F matrices
    counti=1
    do i=1,d
        do j=1,d
            if (i /= j) then
                b1%matrix=f(i,j,d)
                call writenum(i,c,'i2')
                call writenum(j,c2,'i2')
                call removesp(c)
                call removesp(c2)
                b1%idx='f('//c//','//c2//')'
                call removesp(b1%idx)
                rrbas%bas(d+counti)=b1
                counti=counti+1
            end if
        end do
    end do
    rrbas%dimen=d
    print *,'Done generating basis.'
end function genHbasis

function genUbasis(d) result(rrbas)
    integer, intent(in) :: d
    type(Ubas) :: b1
    type(Ubasis) :: rrbas
    integer :: i,j,counti,ed
    character(len=2) c,c1,c2

    ed=d**2
    allocate(rrbas%bas(ed))

    print *,'Generating unitary basis ...'

    do i=1,ed
        allocate(rrbas%bas(i)%matrix(d,d))
    end do

    call writenum(d,c1,'i2')
    call removesp(c1)
    counti=1

    do i=1,d
        do j=1,d
                b1%matrix=genPauli(i,j,d)
                call writenum(i,c,'i2')
                call writenum(j,c2,'i2')
                call removesp(c)
                call removesp(c2)
                b1%idx='S('//c//','//c2//')'
                call removesp(b1%idx)
                rrbas%bas(counti)=b1
                counti=counti+1
        end do
    end do
    rrbas%dimen=d
    print *,'Done generating basis.'
end function genUbasis

subroutine printHbasis(bb)
    type(Hbasis), intent(in) :: bb
    integer i,ed

    ed=perm(bb%dimen,2)+bb%dimen

    do i=1,ed
        call DISP(bb%bas(i)%idx//' = ',bb%bas(i)%matrix)
    end do

    return
end subroutine printHbasis

subroutine printUbasis(bb)
    type(Ubasis), intent(in) :: bb
    integer i,ed

    ed=bb%dimen**2

    do i=1,ed
        call DISP(bb%bas(i)%idx//' = ',bb%bas(i)%matrix)
    end do

    return
end subroutine printUbasis

function cart3d_to_h2(x,y,z) result(rr)
    type(components) :: rr
    real, intent(in) :: x,y,z
    real :: nrm

    nrm=sqrt(x*x+y*y+z*z)

    rr%cmp(1)%idx='f(1,2)'
    rr%cmp(1)%x=x/nrm
    rr%cmp(2)%idx='f(2,1)'
    rr%cmp(2)%x=y/nrm
    rr%cmp(3)%idx='f(2,2)'
    rr%cmp(3)%x=z/nrm

end function cart3d_to_h2


end module