program nbody
 use config
 use step
 use force
 implicit none
 real, dimension(ndim) :: x,v,f,x0,v0
 real :: dt,e,t,dtsub,tsub,deltat
 integer :: i,j,k,nsteps,nf,nsub
 integer :: iunit(nmethods), ievunit(nmethods)
 logical :: active
 
 do j=1,nmethods
    open(newunit=iunit(j),  file=trim(method(j))//'.out',status='replace')
    open(newunit=ievunit(j),file=trim(method(j))//'.ev',status='replace')
 enddo

 do j=1,nmethods
    e = 0.7
    x = (/1.-e,0./)
    v = (/0.,sqrt((1.+e)/(1.-e))/)
    x0 = x
    v0 = v
    t = 0.
    dt = 0.1
    nsteps = nint(200./dt)
    nf = 0 ! number of force evaluations
    
    call get_force(f,x)
    write(iunit(j),*) x
    write(ievunit(j),*) t,nf,x(1)*v(2) - x(2)*v(1),0.5*dot_product(v,v) - 1./norm2(x)

    do i=1,nsteps
       active = .true.
       select case(trim(method(j)))
       case('pefrl')
          call step_pefrl(x,v,f,dt,nf)
       case('fr4')
          call step_fr4(x,v,f,dt,nf)
       case('verlet')
          call step_verlet(x,v,f,dt,nf,active)
       case('verletind')
          nsub = 8
          dtsub = dt/nsub
          tsub  = t
          do k=1,nsub
             tsub = tsub + dtsub
             active = (k==8)
             deltat = tsub - t
             !print*, 'dtsub = ',k,deltat,dt,deltat/dt
             call step_verlet(x,v,f,dtsub,nf,active)
             !call step_verletind(x,v,f,deltat,nf,active,x0,v0)
          enddo
       case('rk2')
          call step_rk2(x,v,f,dt,nf)
       case('rk4')
          call step_rk4(x,v,f,dt,nf)
       case('leapfrog')
          call step_leapfrog(x,v,f,dt,nf)
       end select
       t = t + dt
       write(iunit(j),*)   x
       write(ievunit(j),*) t,nf,x(1)*v(2) - x(2)*v(1),0.5*dot_product(v,v) - 1./norm2(x)
    enddo
 enddo

 do j=1,nmethods
    close(iunit(j))
    close(ievunit(j))
 enddo

end program nbody
