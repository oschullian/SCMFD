      !
      !  Purpose: 
      !    This subroutine generates the velocity of the field particle
      !    with whom the particle of interest is colliding
      !    
      !    
      !
      !
      !  Date       Programmer        Description of change
      !  ====       ==========        =====================
      !  04/05/14   Otto Schullian    
      !


      subroutine generatecolvel(pt, v2, ind)
      use param
      implicit none
      
      ! the type of the field 
      integer, intent(in) :: pt
      ! output velocity of the field particle
      double precision, dimension(3), intent(out) :: v2
      ! index in the avp array of the velocity v2
      integer, intent(out) :: ind
      integer :: ix, iy, iz, ptype
      
      
      ! random number for random (when mkl arrives, change this)
      ! real, dimension(nopco) :: rvec
      ! loop variable
      integer :: k
      ! relative velocities of all particles of the field
      double precision, dimension(noavp) :: vrel
      ! random value to decide which field particle is chosen
      real, dimension(1) :: rvecrandvalue
      double precision, dimension(1) :: randvalue
      
      double precision :: ccs, r1, r2
      
      ! retrieve subcell data and particle type 
      ix=nint(particle(8))
      iy=nint(particle(9))
      iz=nint(particle(10))
      ptype=nint(particle(7))
      ! initialise vrel
      vrel=0.d0
      
      
      ! for all these random positions calculate the the relative velocity 
      ! with respect to vr, immediately sum up to obtain cdf
      vrel(1)=sum((avpvel3(1,pt,:,ix,iy,iz)-particle(4:6))**2)
      
      if (vrel(1)==0.d0) then
        vrel(1)=0.d0
      else
      ! calculate radius of both particle types
      r2=ccsfa(ptype,pt)/
     &dsqrt(dexp((visind(pt)-0.5d0)*dlog(vrel(1))))/2.d0
      r1=ccsfa(pt,ptype)/
     &dsqrt(dexp((visind(ptype)-0.5d0)*dlog(vrel(1))))/2.d0
      ! calculate cross section
      ccs=pi*(r1+r2)**2
      ! calculate vrel
      vrel(1)=ccs*dsqrt(vrel(1))
      end if
      
      ! do this for all particles 
      do k=2,nopcopos3(pt,ix,iy,iz)
        vrel(k)=
     &   sum((avpvel3(k,pt,:,ix,iy,iz)-particle(4:6))**2)

        if (vrel(k)==0.d0) then
          vrel(k)=vrel(k-1)
        else
        
        r2=ccsfa(ptype,pt)/
     &dsqrt(dexp((visind(pt)-0.5d0)*dlog(vrel(k))))/2.d0
        r1=ccsfa(pt,ptype)/
     &dsqrt(dexp((visind(ptype)-0.5d0)*dlog(vrel(k))))/2.d0
        ccs=pi*(r1+r2)**2

          vrel(k)=vrel(k-1)+ccs*dsqrt(vrel(k))
        end if
      end do
      
      ! if all particles have the same velocity then the vector vrel is zero
      ! theoretically no collision occurs, but since vr is zero, the collision
      ! routine will not change the velocities, this means that we can return the first
      ! entry of the vector
      if (vrel(nopcopos3(pt,ix,iy,iz)).le.0.e0) then
       v2=avpvel3(1,pt,:,ix,iy,iz)
       ind=1
       return
      end if 
      
      ! normalise the cdf
      vrel=vrel/vrel(nopcopos3(pt,ix,iy,iz))

      
      
      ! generate one number to decide which particle collides
      call RANMAR(rvecrandvalue, 1)
      randvalue=dble(rvecrandvalue)
     
      ! decide which particle collides
      do k=1,nopcopos3(pt,ix,iy,iz)
        if (randvalue(1).le.vrel(k)) then
          ! save the index of the particle
          ind=k
          exit
        end if
      end do
      
      
     
      
      ! save the velocity of the field particle
      v2=avpvel3(ind,pt,:,ix,iy,iz)

      
      
      end subroutine generatecolvel
      
