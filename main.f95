! main prog
program main
use basis
use dispmodule
use fact
implicit none

type(Hbasis) :: bb1
type(Ubasis) :: bb2

call DISP('f = ',f(3,2,3))

call DISP('dirsum = ',dirsum(f(2,3,3),3))

allocate(bb1%bas(perm(4,2)+4))
bb1=genHbasis(4)

call printBasis(bb1)

allocate(bb2%bas(4**2))
bb2=genUbasis(4)

call printBasis(bb2)

end program