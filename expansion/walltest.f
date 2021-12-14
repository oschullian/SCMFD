      !
      !  Purpose: 
      !    This subroutine test if a particle has crashed into a wall when crossing 
      !    the boundary of a subcell. If not it returns otherwise the wallhit routine
      !    is called
      !    
      !    
      !
      !
      !  Date       Programmer        Description of change
      !  ====       ==========        =====================
      !  18/11/19   Otto Schullian    
      !


      subroutine walltest(index)
      use param
      implicit none
      
      ! the dimension through which the particle has crossed a subcell boundary
      integer, intent(in) :: index
      ! the cell the particle was in
      integer :: pcell
      ! testing parameter if particle is still within the cell
      integer :: testval
      
      
      ! determine the particle subcell
      pcell=nint(particle(7+index))
      
      ! calculate this to check if particle has a subcell 0 or subcell ncell+1
      testval=pcell*(pcell-ncellall(index)-1)
      
      ! if it is inside, return
      if (testval<0.d0) return
      

      ! if it is the wall in z-direction, just loop it back
      if (index==1) then
        if (pcell==0) then
          call hitwall(index)
          
        else
        ! if the wall is in y or z-direction define the trajectory as over and return
        parout=1
       
        end if
      else
        ! if the wall is in y or z-direction define the trajectory as over and return
        parout=1
      
      end if
      
      
      
      end subroutine walltest
      

      
      
