module force
 use config
 implicit none

contains

 subroutine get_force(f,x)
  real, dimension(ndim), intent(in)  :: x
  real, dimension(ndim), intent(out) :: f
  real :: r,r2
  
  r2 = dot_product(x,x)
  r = sqrt(r2)
  f = -x/(r2*r)

 end subroutine get_force

end module force
