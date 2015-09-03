! basis creation
!
module basis
use types
use constants, only: I
use strings
use fact
use dispmodule
implicit none

type Hbas
    complex(dp), dimension(:,:), allocatable :: matrix
    character(8) idx
    logical :: hermitian = .true.
    logical :: unitary = .false.
end type Hbas

type Ubas
    complex(dp), dimension(:,:), allocatable :: matrix
    character(7) idx
    logical :: unitary = .true.
    logical :: hermitian = .false.
end type Ubas

type HBasis
    type(Hbas), dimension(:), allocatable :: bas
    integer dimen
end type HBasis

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

subroutine printHbasis(bb)
    type(Hbasis), intent(in) :: bb
    integer i,ed

    ed=perm(bb%dimen,2)+bb%dimen

    do i=1,ed
        call DISP(bb%bas(i)%idx//' = ',bb%bas(i)%matrix)
    end do

    return
end subroutine
end module