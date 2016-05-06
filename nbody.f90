module utils
 implicit none
 integer, parameter :: ndim = 2

contains

 subroutine get_force(f,x)
  real, dimension(ndim), intent(in)  :: x
  real, dimension(ndim), intent(out) :: f
  real :: r,r2
  
  r2 = dot_product(x,x)
  r = sqrt(r2)
  f = -x/(r2*r)

 end subroutine get_force
 
 subroutine step_verlet(x,v,f,dt)
  real, dimension(ndim), intent(inout) :: x,v,f
  real, intent(in) :: dt
  real, dimension(ndim) :: fprev

  fprev = f
  x = x + dt*v + 0.5*dt**2*f
  call get_force(f,x)
  v = v + 0.5*dt*(f + fprev)

 end subroutine step_verlet

 subroutine step_leapfrog(x,v,f,dt)
  real, dimension(ndim), intent(inout) :: x,v,f
  real, intent(in) :: dt

  x = x + 0.5*dt*v
  call get_force(f,x)
  v = v + dt*f
  x = x + 0.5*dt*v

 end subroutine step_leapfrog

 subroutine step_rk2(x,v,f,dt)
  real, dimension(ndim), intent(inout) :: x,v,f
  real, intent(in) :: dt
  real, dimension(ndim) :: xin,vin

  xin = x
  vin = v
  x = xin + 0.5*dt*v
  v = vin + 0.5*dt*f
  call get_force(f,x)
  x = xin + dt*v
  v = vin + dt*f
  call get_force(f,x)

 end subroutine step_rk2

 subroutine step_rk4(x,v,f,dt)
  real, dimension(ndim), intent(inout) :: x,v,f
  real, intent(in) :: dt
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

 end subroutine step_rk4

 !--Forest-Ruth 4th order symplectic
 subroutine step_fr4(x,v,f,dt)
  real, dimension(ndim), intent(inout) :: x,v,f
  real, intent(in) :: dt
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

 end subroutine step_fr4

 !--Position-Extended Forest-Ruth Like (PEFRL) 4th order symplectic
 !  (no backwards steps, but 4 force evaluations)
 subroutine step_pefrl(x,v,f,dt)
  real, dimension(ndim), intent(inout) :: x,v,f
  real, intent(in) :: dt
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

 end subroutine step_pefrl

end module utils

program kepler
 use utils
 implicit none
 real, dimension(ndim) :: x,v,f
 real :: dt,e,t
 integer :: i,j,nsteps
 integer, parameter :: ntry = 5
 
 open(unit=1,file='lf.out',status='replace')
 open(unit=2,file='rk.out',status='replace')
 open(unit=3,file='verlet.out',status='replace')
 open(unit=4,file='fr4.out',status='replace')
 open(unit=5,file='pefrl.out',status='replace')
 open(unit=6,file='lf.ev',status='replace')
 open(unit=7,file='rk.ev',status='replace')
 open(unit=8,file='verlet.ev',status='replace')
 open(unit=9,file='fr4.ev',status='replace')
 open(unit=10,file='pefrl.ev',status='replace')

 do j=1,ntry
    e = 0.7
    x = (/1.-e,0./)
    v = (/0.,sqrt((1.+e)/(1.-e))/)
    t = 0.
    dt = 0.1
   ! nsteps = 100./dt
    
    call get_force(f,x)
    write(j,*) x
    write(ntry+j,*) t,x(1)*v(2) - x(2)*v(1),0.5*dot_product(v,v) - 1./norm2(x)

    !if (j.eq.1) x = x + 0.25*f*dt**2

    do i=1,5000
       select case(j)
       case(5)
          call step_pefrl(x,v,f,dt)       
       case(4)
          call step_fr4(x,v,f,dt)       
       case(3)
          call step_verlet(x,v,f,dt)       
       case(2)
          call step_rk4(x,v,f,dt)
       case(1)
          call step_leapfrog(x,v,f,dt)
       end select
       t = t + dt
       write(j,*) x
       write(ntry+j,*) t,x(1)*v(2) - x(2)*v(1),0.5*dot_product(v,v) - 1./norm2(x)
    enddo
 enddo

 do j=1,2*ntry
    close(j)
 enddo

end program kepler
