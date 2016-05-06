program nbody
 use config
 use step
 use force
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
    nsteps = 5000 !100./dt
    
    call get_force(f,x)
    write(j,*) x
    write(ntry+j,*) t,x(1)*v(2) - x(2)*v(1),0.5*dot_product(v,v) - 1./norm2(x)

    !if (j.eq.1) x = x + 0.25*f*dt**2

    do i=1,nsteps
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

end program nbody
