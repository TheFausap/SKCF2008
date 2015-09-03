! main prog
program main
use basis
use dispmodule
use fact
implicit none

type(Hbasis) :: bb1
call DISP('f = ',f(3,2,3))

call DISP('dirsum = ',dirsum(f(2,3,3),3))

!call DISP('h(3,3) = ',h(3,3))

allocate(bb1%bas(perm(4,2)+4))
bb1=genHbasis(4)
!call DISP(bb1%bas(3)%idx//' = ',bb1%bas(3)%matrix)

!call DISP(bb1%bas(2)%idx//' = ',bb1%bas(2)%matrix)

call printHbasis(bb1)

end program