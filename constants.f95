module constants
use types, only: dp
implicit none
private
public pi, e, I
! Constants contain more digits than double precision, so that
! they are rounded correctly:
real(dp), parameter :: pi   = 3.1415926535897932384626433832795_dp
real(dp), parameter :: e    = 2.7182818284590452353602874713527_dp
complex(dp), parameter :: I = (0, 1)
end module