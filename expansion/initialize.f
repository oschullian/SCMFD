      !
      !  Purpose: 
      !    This subroutine initializes all the initial parameters
      !    
      !    
      !
      !
      !  Date       Programmer        Description of change
      !  ====       ==========        =====================
      !  04/05/14   Otto Schullian    
      !


      subroutine initialize
      use param
      implicit none
      
      
      character(len=100) :: dum
      
      integer :: readerror
      integer :: k,l
      
      double precision :: val, mured, ccs, dtmin, vel, redm, pres, dens
      double precision, dimension(:), allocatable :: dtprinttest
      double precision, dimension(:), allocatable :: ui
      !initialise value of pi
      pi=acos(-1.e0) 
      
      
      !open the input file and check whether opening was successful
      open(unit=1, file='input.txt', action='READ', status='old',
     & iostat=readerror)
      if (readerror.ne.0) then 
        write(*,*) 'can''t open source file'
        write(*,*) 'check path'
        stop
      end if
      
      !discard the first line
      read(1,*,iostat=readerror) dum
      
      write(*,*) dum
      
      !number of particle types
      read(1,*,iostat=readerror) dum(1:3), nop
      if (readerror.ne.0) then 
        write(*,*) 'parameter read error'
        stop
      end if
      
      write(*,*) 'nop', nop
      
      ! number of particles for averages
      read(1,*,iostat=readerror) dum(1:5), noavp
      if (readerror.ne.0) then 
        write(*,*) 'parameter read error'
        stop
      end if
      
      write(*,*) 'noavp', noavp
      
      
 
      
      ! number of total runs
      read(1,*,iostat=readerror) dum(1:7), totruns
      if (readerror.ne.0) then 
        write(*,*) 'parameter read error'
        stop
      end if
      
      write(*,*) 'totruns', totruns
      
      
      

      ! avpstart
      read(1,*,iostat=readerror) dum(1:8), avpstart
      if (readerror.ne.0) then 
        write(*,*) 'parameter read error'
        stop
      end if
      
      write(*,*) 'avpstart', avpstart
      nprint=avpstart
      nprintdens=avpstart
      
      ! modavpout
      read(1,*,iostat=readerror) dum(1:6), avpout
      if (readerror.ne.0) then
        write(*,*) 'parameter read error'
        stop
      end if
      
      write(*,*) 'avpout', avpout

      ! modavpout
      read(1,*,iostat=readerror) dum(1:9), modmomout
      if (readerror.ne.0) then
        write(*,*) 'parameter read error'
        stop
      end if
      
      write(*,*) 'modmomout', modmomout
      ! modavpout
      read(1,*,iostat=readerror) dum(1:9), modnewavp
      if (readerror.ne.0) then
        write(*,*) 'parameter read error'
        stop
      end if
      
      write(*,*) 'modnewavp', modnewavp
      
      ! modavpout
      read(1,*,iostat=readerror) dum(1:11), inpcellperc
      if (readerror.ne.0) then 
        write(*,*) 'parameter read error'
        stop
      end if
      
      write(*,*) 'inpcellperc', inpcellperc
      
      
      
      
      ! modavpout
      read(1,*,iostat=readerror) dum(1:4), seed
      if (readerror.ne.0) then 
        write(*,*) 'parameter read error'
        stop
      end if
      
      write(*,*) 'seed', seed

      ! modavpout
      read(1,*,iostat=readerror) dum(1:10), fileappend
      if (readerror.ne.0) then
        write(*,*) 'parameter read error'
        stop
      end if

      write(*,*) 'fileappend', fileappend
      
      ! modavpout
      read(1,*,iostat=readerror) dum(1:10), radoutflow
      if (readerror.ne.0) then
      write(*,*) 'parameter read error'
      stop
      end if

      write(*,*) 'radoutflow', radoutflow
   
      
      ! dummy
      read(1,*,iostat=readerror) dum
      if (readerror.ne.0) then 
        write(*,*) 'parameter read error'
        stop
      end if
      
      ! distance between walls
      read(1,*,iostat=readerror) dum(1:5), distx
      if (readerror.ne.0) then 
        write(*,*) 'parameter read error'
        stop
      end if
      read(1,*,iostat=readerror) dum(1:5), disty
      if (readerror.ne.0) then 
        write(*,*) 'parameter read error'
        stop
      end if
      read(1,*,iostat=readerror) dum(1:5), distz
      if (readerror.ne.0) then 
        write(*,*) 'parameter read error'
        stop
      end if
      
      write(*,*) 'distance', distx, disty, distz
      distall(1)=distx
      distall(2)=disty
      distall(3)=distz
      
      ! number of cells between walls
      read(1,*,iostat=readerror) dum(1:5), ncellx
      if (readerror.ne.0) then 
        write(*,*) 'parameter read error'
        stop
      end if
      ! number of cells between walls
      read(1,*,iostat=readerror) dum(1:5), ncelly
      if (readerror.ne.0) then 
        write(*,*) 'parameter read error'
        stop
      end if
      ! number of cells between walls
      read(1,*,iostat=readerror) dum(1:5), ncellz
      if (readerror.ne.0) then 
        write(*,*) 'parameter read error'
        stop
      end if
      ncellall(1)=ncellx
      ncellall(2)=ncelly
      ncellall(3)=ncellz
      
      write(*,*) 'ncell', ncellx, ncelly, ncellz
      
      ! number of cells between walls
      read(1,*,iostat=readerror) dum(1:8), tempwall
      if (readerror.ne.0) then 
        write(*,*) 'parameter read error'
        stop
      end if
      
      write(*,*) 'tempwall', tempwall
      
      
      ! number of cells between walls
      read(1,*,iostat=readerror) dum(1:11), radiusvalve
      if (readerror.ne.0) then 
        write(*,*) 'parameter read error'
        stop
      end if
      if (radiusvalve<0.d0) then
        write(*,*) 'radius is negative'
        stop
      end if 
      write(*,*) 'radiusvalve', radiusvalve
      
      
      
      
      ! allocate arrays that contain information on particles
      allocate(np(nop))
      allocate(mass(nop))
      allocate(influx(nop))
      
      allocate(refdia(nop))
      allocate(reftemp(nop))
      allocate(visind(nop))
      allocate(dtprint(nop))
      ! loop over particle parameters
      do k=1,nop
        write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(*,*) 'parameters particle', k
        read(1,*,iostat=readerror) dum
        if (readerror.ne.0) then 
          write(*,*) 'parameter read error'
          stop
        end if
        ! number of trajectories per run
        read(1,*,iostat=readerror) dum(1:2), np(k)
        if (readerror.ne.0) then 
          write(*,*) 'parameter read error'
          stop
        end if
      
        write(*,*) 'np', np(k)
        
      
      
        ! particle reference diameter
        read(1,*,iostat=readerror) dum(1:4), refdia(k)
        if (readerror.ne.0) then 
          write(*,*) 'parameter read error'
          stop
        end if
      
        write(*,*) 'refdia', refdia(k)
        
    
      
        ! particle reference temperature 
        read(1,*,iostat=readerror) dum(1:4), reftemp(k)
        if (readerror.ne.0) then 
          write(*,*) 'parameter read error'
          stop
        end if
      
        write(*,*) 'reftemp', reftemp(k)
        
        ! viscosity index 
        read(1,*,iostat=readerror) dum(1:4), visind(k)
        if (readerror.ne.0) then 
          write(*,*) 'parameter read error'
          stop
        end if
      
        write(*,*) 'visind', visind(k)
        
        
        ! mass
        read(1,*,iostat=readerror) dum(1:9), mass(k)
        if (readerror.ne.0) then 
          write(*,*) 'parameter read error'
          stop
        end if
      
        write(*,*) 'mass', mass(k)
        
        ! particle density
        read(1,*,iostat=readerror) dum(1:6), influx(k)
        if (readerror.ne.0) then 
          write(*,*) 'parameter read error'
          stop
        end if
      
        write(*,*) 'influx', influx(k)
        ! particle density
          
      end do

      allocate(dtprinttest(nop))
      allocate(ui(nop))
      write(*,*) '----------------------------------------'
      write(*,*) 'checking if dtprint is too high or low'
      write(*,*) 'the average time a particle travels'
      write(*,*) 'in the cell is given by'
      ui=dsqrt(5.d0*AK*tempwall/(mass))
      write(*,*) maxval(distall)
     &  /ui
      write(*,*) 'number of particles flown before avp is updated'
      write(*,*) dble(np)*dble(modnewavp)
      write(*,*) 'the total time of all particles spent in the'
      write(*,*) 'system is'
      write(*,*) maxval(distall)
     &  /ui
     &  *dble(np)*dble(modnewavp)
      write(*,*) 'the number of cells'
      write(*,*) product(ncellall),ncellall
      write(*,*) 'the number of updates that should happen'
      write(*,*) product(ncellall)*inpcellperc*noavp
      write(*,*) 'the number of updates that should happen p.p'
      write(*,*) product(ncellall)*inpcellperc*noavp/dble(np)
      write(*,*) 'to acchieve this n.o.updates dtprint should be'
      dtprinttest=maxval(distall)
     &  /ui
     &  *dble(np)*dble(modnewavp)
     &  /(product(ncellall)*inpcellperc*noavp)
      write(*,*) dtprinttest
      dtprint=dtprinttest/2.d0
      write(*,*) 'dtprint is'
      write(*,*) dtprint
      write(*,*) 'with the current dtprint this n.o updates are'
      write(*,*) 'in total'
      write(*,*) maxval(distall)
     &  /ui
     &  *dble(np)*dble(modnewavp)/dtprint
      write(*,*) 'p.p'
       write(*,*) maxval(distall)
     &  /ui
     &  /dtprint
      write(*,*) '----------------------------------------'
      



 



      close(unit=1)
      
      
      ! allocate arrays that hold simulation information
      ! array that contains velocities of previous particles
      allocate(avpvel3(noavp,nop,3,ncellx,ncelly,ncellz))
      ! allocate(avpveltest(noavp,nop,3,ncell))
      ! counter to find the next position in array avpvel to save
      ! velocity
      allocate(avtexchange3(nop,ncellx,ncelly,ncellz))
      
      allocate(avtpos3(nop,ncellx,ncelly,ncellz))
     
      !allocate(avpvelpos(nop,ncell))
      !allocate(avpvelposprev(nop,ncell))
      !allocate(avpvelposnew(nop,ncell))
      allocate(nopcopos3(nop,ncellx,ncelly,ncellz))

      !allocate(nopcoposnew(nop,ncell))
      ! array for calculation of moments
      allocate(n_count3(nop,ncellx,ncelly,ncellz))
      ! array for calculation of first moment
      allocate(mom_fi3(3,nop,ncellx,ncelly,ncellz))
      ! array for calculation of second moment
      allocate(mom_sec3(6,nop,ncellx,ncelly,ncellz))
      allocate(n_countav(nop,ncellx,ncelly,ncellz))
      ! array for calculation of first moment
      allocate(mom_fiav(3,nop,ncellx,ncelly,ncellz))
      ! array for calculation of second moment
      allocate(mom_secav(6,nop,ncellx,ncelly,ncellz))
      ! kval of previous calculations
      allocate(kvalprev(nop))
      ! vector containing information on kval calculation
      ! allocate(kvalcalc(nop))
      n_count3=0.d0
      mom_fi3=0.d0
      mom_sec3=0.d0
      n_countav=0.d0
      mom_fiav=0.d0
      mom_secav=0.d0
    
      allocate(ncelldens3(nop,ncellx,ncelly,ncellz))
      allocate(ncelldensnew3(nop,ncellx,ncelly,ncellz))
     
      
      
      ! calculate general parameters
      allocate(ccsfa(nop,nop))
      do k=1,nop
        
        call gammacalc(visind(k),val)
        
        do l=1,nop
          mured=1.d0/mass(k)+1.d0/mass(l)
          mured=1.d0/mured
          
          ccsfa(k,l)=refdia(k)*dsqrt(dexp((visind(k)-0.5d0)*
     &      dlog(2.d0*AK*reftemp(k)/
     &      mured))/val)
          
          write(*,*)  k,l, ccsfa(k,l)
        end do
        
      end do 
      
      write(*,*) 'ccsfa', ccsfa
      
      
      
      
      !levicivita=0.e0
      !levicivita(1,2,3)=1.e0
      !levicivita(2,3,1)=1.e0
      !levicivita(3,1,2)=1.e0
      !levicivita(1,3,2)=-1.e0
      !levicivita(3,2,1)=-1.e0
      !levicivita(2,1,3)=-1.e0
      
      
      
      
      
      allocate(colcount(nop,nop,ncellx,ncelly,ncellz))
      allocate(disttrav(nop))
      allocate(timecounter(nop))
      colcount=0.d0
      disttrav=0.d0
      timecounter=0.d0



      nvelchange=ceiling(inpcellperc*dble(noavp))
      write(*,*) 'nvelchange', nvelchange
      allocate(avpvelinter3(nvelchange,nop,3,ncellx,ncelly,ncellz))



      
      
      write(*,*) '--------------------------------------------'
      write(*,*) 'the traj time in respect to the cellsize is evaluated'
      
      do k=1,nop
          write(*,*) 'particle', k
          write(*,*) 'cellsize', distx/dble(ncellx),
     &    disty/dble(ncelly),distz/dble(ncellz)
          
          vel=dsqrt(8.d0*AK*tempwall/(mass(k)*pi))
          
          write(*,*) 'dt','=', distx/dble(ncellx)/vel,
     & disty/dble(ncelly)/vel, distz/dble(ncellz)/vel
    
      end do
      write(*,*) '--------------------------------------------'
      
      write(*,*) '--------------------------------------------'
      write(*,*) 'the fluxes correspond to the following pressures:'
      do k=1,nop
        write(*,*) 'for particle: ', k
        pres=influx(k)*
     &  dsqrt(2.d0*pi*mass(k)*mconv*AK2*tempwall)
     &  /(pi*radiusvalve**2)
        dens=influx(k)*
     &  dsqrt(2.d0*pi*mass(k)*mconv/(AK2*tempwall))
     &  /(pi*radiusvalve**2)
        write(*,*) 'pressure = ', pres, 'Pa'
        write(*,*) 'pressure = ', pres/1.d5, 'bar'
        write(*,*) 'density = ', dens, 'm-3'
      end do 
      
      
      
      write(*,*) '--------------------------------------------'
      
      
      
      
      
      end subroutine initialize
      

      
      
