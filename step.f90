module step
 use config
 use force
 implicit none
 integer, parameter :: nmethods = 6
 character(len=*), parameter :: method(nmethods) = &
   (/'verlet  ', &
     'leapfrog',&
     'rk2     ',&
     'rk4     ',&
     'fr4     ',&
     'pefrl   '/)

contains
 
 subroutine step_verlet(x,v,f,dt,nf)
  real, dimension(ndim), intent(inout) :: x,v,f
  real, intent(in) :: dt
  integer, intent(inout) :: nf
  real, dimension(ndim) :: fprev

  fprev = f
  x = x + dt*v + 0.5*dt**2*f
  call get_force(f,x)
  v = v + 0.5*dt*(f + fprev)  
  nf = nf + 1

 end subroutine step_verlet

 subroutine step_leapfrog(x,v,f,dt,nf)
  real, dimension(ndim), intent(inout) :: x,v,f
  real, intent(in) :: dt
  integer, intent(inout) :: nf

  x = x + 0.5*dt*v
  call get_force(f,x)
  v = v + dt*f
  x = x + 0.5*dt*v
  nf = nf + 1

 end subroutine step_leapfrog

 subroutine step_rk2(x,v,f,dt,nf)
  real, dimension(ndim), intent(inout) :: x,v,f
  real, intent(in) :: dt
  integer, intent(inout) :: nf
  real, dimension(ndim) :: xin,vin

  xin = x
  vin = v
  x = xin + 0.5*dt*v
  v = vin + 0.5*dt*f
  call get_force(f,x)
  x = xin + dt*v
  v = vin + dt*f
  call get_force(f,x)
  nf = nf + 2

 end subroutine step_rk2

 subroutine step_rk4(x,v,f,dt,nf)
  real, dimension(ndim), intent(inout) :: x,v,f
  real, intent(in) :: dt
  integer, intent(inout) :: nf
  real, dimension(ndim) :: xin,vin
  real, dimension(ndim) :: v2,v3,v4,f2,f3,f4

  xin = x
  vin = v

  x = xin + 0.5*dt*v
  v = vin + 0.5*dt*f

  v2 = v
  call get_force(f2,x)

  x = xin + 0.5*dt*v2
  v = vin + 0.5*dt*f2

  v3 = v
  call get_force(f3,x)

  x = xin + dt*v3
  v = vin + dt*f3

  v4 = v
  call get_force(f4,x)  

  x = xin + dt/6.*(vin + 2.*v2 + 2.*v3 + v4)
  v = vin + dt/6.*(f + 2.*f2 + 2.*f3 + f4)

  call get_force(f,x)

  nf = nf + 4

 end subroutine step_rk4

 !--Forest-Ruth 4th order symplectic
 subroutine step_fr4(x,v,f,dt,nf)
  real, dimension(ndim), intent(inout) :: x,v,f
  real, intent(in) :: dt
  integer, intent(inout) :: nf
  real, parameter :: theta = 1./(2. - 2**(1./3.)) 

  x = x + theta*0.5*dt*v
  call get_force(f,x)
  v = v + theta*dt*f
  x = x + (1. - theta)*0.5*dt*v
  call get_force(f,x)
  v = v + (1. - 2.*theta)*dt*f   ! note this is a backwards step...
  x = x + (1. - theta)*0.5*dt*v
  call get_force(f,x)
  v = v + theta*dt*f
  x = x + theta*0.5*dt*v

  nf = nf + 3

 end subroutine step_fr4

 !--Position-Extended Forest-Ruth Like (PEFRL) 4th order symplectic
 !  (no backwards steps, but 4 force evaluations)
 subroutine step_pefrl(x,v,f,dt,nf)
  real, dimension(ndim), intent(inout) :: x,v,f
  real, intent(in) :: dt
  integer, intent(inout) :: nf
  real, parameter :: eps = 0.1786178958448091
  real, parameter :: lam = -0.2123418310626054
  real, parameter :: chi = -0.06626458266981849

  x = x + eps*dt*v
  call get_force(f,x)
  v = v + (1. - 2.*lam)*0.5*dt*f
  x = x + chi*dt*v
  call get_force(f,x)
  v = v + lam*dt*f
  x = x + (1. - 2.*(chi + eps))*dt*v
  call get_force(f,x)
  v = v + lam*dt*f
  x = x + chi*dt*v
  call get_force(f,x)
  v = v + (1. - 2.*lam)*0.5*dt*f
  x = x + eps*dt*v

  nf = nf + 4

 end subroutine step_pefrl

end module step
