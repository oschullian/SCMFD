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


      subroutine hitwall(index)
      use param
      implicit none
      

      integer, intent(in) :: index
      integer :: walltype
      integer ::  pt, pcell
      
      double precision, dimension(3) :: vel
      
      ! retriebe particle type and the cell, that is too high/or low
      pt=nint(particle(7))
      pcell=nint(particle(7+index))
      
      ! if pcell is 0, it will go to 1 and otherwise -1
      if (pcell==0) then
        walltype=1
        particle(7+index)=particle(7+index)+1.d0
      else 
        walltype=2
        particle(7+index)=particle(7+index)-1.d0
      end if

      ! update the wallhit counter
      wallhitstat(index,walltype)=wallhitstat(index,walltype)+1
      
      
      
      ! generate a new wall velocity 
      call generatewallvel(pt,walltype,vel,index)
      
      ! save the new velocity
      particle(4:6)=vel
      
      
      
300   format(3000E23.13)
      
      end subroutine hitwall
      

      
      
