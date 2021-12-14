      !
      !  Purpose: 
      !    This subroutine updates the avpvel whenever called
      !    
      !    
      !
      !
      !  Date       Programmer        Description of change
      !  ====       ==========        =====================
      !  04/05/14   Otto Schullian    
      !


      subroutine timestepprop(time)
      use param
      implicit none
      
      ! the time corresponding to the 
      double precision, intent(inout) :: time
      ! sum for the calculation of the free time
      double precision :: dummysum
      ! time(s) for the exiting the subcell and the free time
      double precision, dimension(5) :: deltat
      
      ! vector for uniform random number
      real, dimension(1) :: rvecwall
      ! the time the particle is actually propagated
      double precision :: tfinal
      
      ! loop index
      integer :: k
      ! particle integer, particle loop integer and indexes for saving the velocity for update
      integer :: pt, pt2, id1, id2
      ! subcell index and the index of deltat for either exiting the subcell or col
      integer :: ix, iy, iz, tindex
      
      double precision :: a, b, c
      integer :: revcheck
      
      deltat=0.d0
      deltat(5)=dtmax
      ! calculate the expectation value of <v_r sigma>
      call calckval
      
      ! calculate the time until the next collision occurs
      call RANMAR(rvecwall,1)

      ! retrieve subcell index
      ix=nint(particle(8))
      iy=nint(particle(9))
      iz=nint(particle(10))
      
      ! retrieve particle type
      pt=nint(particle(7))
      
      ! for the calculation of the free time
      dummysum=0.d0
      do pt2=1,nop
        dummysum=dummysum+kvalprev(pt2)*ncelldens3(pt2,ix,iy,iz)
      end do 
      
      ! of the sum is 0 no collision will occur, so we set the time very high 
      if (dummysum.le.0.d0) then 
      deltat(4)=1.d5
      else
      deltat(4)=-dlog(1.d0-dble(rvecwall(1)))/dummysum
      end if 
      
      ! calculate the time until the particle leaves the subcell
      do k=1,2
      
      ! depending on the sign of the velocity it is either the left (negative) or 
      ! the right plane, if the velocity is zero in that direction we set the time very 
      ! high
      if (particle(3+k)>0.d0) then 
        deltat(k)=(particle(7+k)*distall(k)/dble(ncellall(k))
     &             -particle(k))/particle(3+k)
      elseif (particle(3+k)<0.d0) then 
        deltat(k)=((particle(7+k)-1.d0)*distall(k)/dble(ncellall(k))
     &             -particle(k))/particle(3+k)
      else 
      deltat(k)=1.d5
      end if
      
      ! if delta(k) is negative either something wrong has happened
      if (deltat(k)<0.d0) then
        write(*,*) 'deltat(k) is negative'
        write(*,*) deltat
        write(*,*) particle
        stop
      end if
      
      end do
      
      ! we do the z-direction separate because of the second power in time
      revcheck=0
      ! in case the velocity is positive then it must go through the right layer
      if (particle(6)>0.d0) then
        ! solve with the midnight formula 
        a=0.5d0*accel
        b=particle(6)
        c=particle(3)-particle(10)*distall(3)/dble(ncellall(3))
        deltat(3)=(-b+dsqrt(b**2-4.d0*a*c))/(2.d0*a)
        
      elseif (particle(6)<0.d0) then 
      ! if it is negative it can either reverse the velocity before it crosses the 
      ! wall or not, we must check 
        a=0.5d0*accel
        b=particle(6)
        c=particle(3)-(particle(10)-1.d0)*distall(3)/dble(ncellall(3))
        
        
        if (b**2-4.d0*a*c<0.d0) then 
          revcheck=1
          a=0.5d0*accel
          b=particle(6)
          c=particle(3)-particle(10)*distall(3)/dble(ncellall(3))
          deltat(3)=(-b+dsqrt(b**2-4.d0*a*c))/(2.d0*a)
        else 
          deltat(3)=(-b-dsqrt(b**2-4.d0*a*c))/(2.d0*a)
        end if
      
      else 
        deltat(k)=1.d5
      end if
      
      if (deltat(3)<0.d0) then
        write(*,*) 'deltat(k) is negative'
        write(*,*) deltat
        write(*,*) particle
        stop
      end if
      
      
      
      ! find the minimum in deltat, that decides what is going to happen first
      tindex=minloc(deltat,1)
      ! the minimum is the time time that will be used for propagation
      tfinal=deltat(tindex)
      
      ! update travelled distance
      if (switchCounter==1) then
       disttrav(pt)=disttrav(pt)+dsqrt(sum(particle(4:6)**2))*tfinal
       timecounter(pt)=timecounter(pt)+tfinal
      end if
      ! update density and moments
      call updatedens(tfinal)
      call updatemom(tfinal)
      
      ! the next index for saving vel
      id1=ceiling((time-dtprintoffset)/dtprint(pt))
      ! find the last index for saving vel
      id2=ceiling((time+tfinal-dtprintoffset)/dtprint(pt))-1
      
      ! if this loop is activated, then time < tprint < time+tfinal, i.e. a "saving
      ! point" was stored
      do k=id1, id2
        call updatevel
      end do 
      
      ! adjust the time
      time=time+tfinal
      
      !propagate the particle
      do k=1,2
        particle(k)=particle(k)+tfinal*particle(3+k)
      end do
      
      particle(3)=particle(3)+tfinal*particle(6)+
     &   0.5d0*accel*tfinal**2
      
      particle(6)=particle(6)+accel*tfinal
      ! if tindex was 4 a collision has occured
      if (tindex==4) then
        ! cal collision routine
        call col
        !write(*,*) 'col'
        !write(*,*) particle
      elseif (tindex==5) then
        ! nothing apart from the propagation happens
        !write(*,*) 'prop'
        !write(*,*) particle
      else 
        
        ! particle has exited the subcell and moved to the next one (it is however 
        ! located at the boundary), so adjust the subcell index 
        if (particle(3+tindex)>0.d0) then 
          particle(7+tindex)=particle(7+tindex)+1.d0
        elseif (particle(3+tindex)<0.d0) then 
          particle(7+tindex)=particle(7+tindex)-1.d0
        else 
          write(*,*) 'velocity is zero but still it crosses the wall?'
          stop
        end if


        ! test if particle is still inside by calling the walltest routine
        call walltest(tindex)
        !write(*,*) 'wall'
        !write(*,*) particle
      end if 
      
      if (isnan(particle(4))) then
        write(*,*) particle
        stop
      end if
      
      
      

      
      end subroutine timestepprop
      

      
      
