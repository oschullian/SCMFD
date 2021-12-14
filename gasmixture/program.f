      program dsmc
      use param
      implicit none
      
      
      ! dummy variable for the loops 
      integer :: pt 
      integer :: nparticle
      integer :: run
      integer :: readerror
      double precision :: time
      
      
      write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      write(*,*) '!'
      write(*,*) '!  DSMCcon Simulation'
      write(*,*) '!'
      write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      
     
      ! generate the files for the output (except avp)
      open(unit=106, file='wallhitstat.txt',
     & iostat=readerror)
      if (readerror.ne.0) then 
        write(*,*) 'can''t open moments file'
        write(*,*) 'check path'
        stop
      end if
      
      open(unit=113, file='testcolenergy.txt',
     & iostat=readerror)
      if (readerror.ne.0) then
      write(*,*) 'can''t open moments file'
      write(*,*) 'check path'
      stop
      end if
      
       open(unit=101, file='disttrav.txt',
     & iostat=readerror)
      if (readerror.ne.0) then 
        write(*,*) 'can''t open disttrav file'
        write(*,*) 'check path'
        stop
      end if
      open(unit=102, file='colcount.txt',
     & iostat=readerror)
      if (readerror.ne.0) then 
        write(*,*) 'can''t open colcount file'
        write(*,*) 'check path'
        stop
      end if
      open(unit=103, file='timecounter.txt',
     & iostat=readerror)
      if (readerror.ne.0) then 
        write(*,*) 'can''t open timecounter file'
        write(*,*) 'check path'
        stop
      end if
      open(unit=104, file='moments.txt',
     & iostat=readerror)
      if (readerror.ne.0) then 
        write(*,*) 'can''t open moments file'
        write(*,*) 'check path'
        stop
      end if

      write(*,*) 'initialize'
      ! call routine that reads in input file
      call initialize
      
      ! initialize seed
      CALL RMARIN(seed,n1random,n2random)
      
      
      write(*,*) 'cdfcalc'
      ! call routine that generates the cdf 
      call calccdf
      
      ! call routine that either generates or reads in initial distributions of 
      ! avp and dens
      call readavp
      call readdens
      
      
      particle=0.d0
      write(*,*) 'start running simulation'
      ! loop over runs, this number should ideally be as large as possible 
      do run=1,totruns
      ! this is the loop over the particle type
      particletype: do pt=1, nop
      
      ! every particle can be run with different accuracy, this is the loop for that
      nparticleloop: do nparticle=1,np(pt)
     
     
      wallhitstat=0
      testenergy=0.d0
      ! generate a particle 
      call generateparticle(pt)
      
      
      time=0.d0
      do while (time<trajlength)
      ! propagate the particle for a timeinterval
      call timestepprop(time)
      end do
      
      
      
      write(106,200) pt, wallhitstat
      write(113,305) dble(pt), testenergy
      end do nparticleloop
      end do particletype
      
      
      write(101,300) disttrav
      write(102,300) colcount
      write(103,300) timecounter
      disttrav=0.d0
      colcount=0.d0
      timecounter=0.d0
      
      write(104,300) n_count3, mom_fi3, mom_sec3
      n_count3=0.d0
      mom_fi3=0.d0
      mom_sec3=0.d0
      
      call flush
      
      
      if (mod(run,modavpout)==0) then
      call printdens
      
      call printavp
    
      end if
      
      
      
      end do
      
      ! close units and leave
      close(unit=101)
      close(unit=102)
      close(unit=103)
      close(unit=104)

      write(*,*) 'bye bye'
      
           
100   format(6A20)
200   format(1000000I20)
300   format(3000E23.13)
305   format(1000000E23.13)
306   format(1000000E23.13)
      
      end program dsmc
