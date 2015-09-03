module fact
   private :: facHelper
contains

recursive function facHelper(n, acc) result(returner)
  implicit none
  integer, intent(in) :: n, acc
  integer :: returner
  if (n <= 1) then
    returner = acc
  else
    returner = facHelper(n - 1, n * acc)
  endif
end function facHelper

function factorial(n)
  implicit none
  integer, intent(in) :: n
  integer :: factorial
  factorial = facHelper(n, 1)
end function factorial

function binom(n,k)
    implicit none
    integer, intent(in) :: n,k
    integer :: binom
    binom=factorial(n)/(factorial(n-k)*factorial(k))
end function binom

function perm(n,k)
    implicit none
    integer, intent(in) :: n,k
    integer :: perm
    perm=factorial(n)/factorial(n-k)
end function perm

end module fact
