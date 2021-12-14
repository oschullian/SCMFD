      !
      !  Purpose: 
      !    This subroutine calculates <vrel sigma > of a particle with a background field
      !    
      !    
      !
      !
      !  Date       Programmer        Description of change
      !  ====       ==========        =====================
      !  04/05/14   Otto Schullian    
      !


      subroutine calckval
      use param
      implicit none
      
      
      ! loop integer
      integer :: k
      ! it is the sum over vrel *sigma, needed for the average
      double precision :: velsum, vrel, r1, r2, ccs, kval
      
      ! particle type of the particle that collides with the field 
      integer :: ptype
      ! particle type of the field (loop index)
      integer ::  pt2
      ! subcellindexes
      integer :: ix, iy, iz

      
      
      ! if the value has to be calculated retrieve type of the particle which 
      ! collides with the field 
      ! retrieve particle type
      ptype=nint(particle(7))
      ix=nint(particle(8))
      iy=nint(particle(9))
      iz=nint(particle(10))
      
      ! loop over all particles
      do pt2=1,nop
      ! initialise value for the sum
      velsum=0.d0
      
      ! loop over all entries in avpvel that have a velocity stored
      do k=1,nopcopos3(pt2,ix,iy,iz)
      
        ! calculate relative velocity 
        vrel=sum((avpvel3(k,pt2,:,ix,iy,iz)-particle(4:6))**2)
        ! if it is zero, then nothing happens
        if (vrel==0.d0) cycle
        
        ! calculate the radii of both particles
        r2=ccsfa(ptype,pt2)/
     &dsqrt(dexp((visind(pt2)-0.5d0)*dlog(vrel)))/2.d0
        r1=ccsfa(pt2,ptype)/
     &dsqrt(dexp((visind(ptype)-0.5d0)*dlog(vrel)))/2.d0
        ! total cross section between these particles
        ccs=pi*(r1+r2)**2
        
        
        ! add vr*ccs to the velsum
          velsum=velsum+
     &   dsqrt(vrel)*ccs
     
      end do
      
      ! divide sum by number of summands in order to obtain the average relative
      ! velocity*cross section
      if (nopcopos3(pt2,ix,iy,iz)==0) then 
      kval=0.d0
      else
      kval=velsum/dble(nopcopos3(pt2,ix, iy, iz))
      end if
      
      ! save the kval in the vector that contains it for all the fields
      kvalprev(pt2)=kval
      end do
      
      
      end subroutine calckval
      
      
      
