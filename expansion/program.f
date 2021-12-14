      program dsmc
      use param
      implicit none
      
      
      ! dummy variable for the loops 
      integer :: pt , ip
      integer :: nparticle
      integer :: run
      integer :: readerror
      double precision :: time
      double precision :: volsubcell, t0, tprev, tnow, distsq
      double precision, dimension(:), allocatable :: runcounter
      
      
      write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      write(*,*) '!'
      write(*,*) '!  DSMCcon Simulation'
      write(*,*) '!'
      write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      
     
      CALL CPU_TIME(t0)
      tprev=t0

      write(*,*) 'initialize'
      ! call routine that reads in input file
      call initialize
      allocate(runcounter(nop))
      runcounter=0.d0
      fileopen: if (fileappend==0) then
      ! generate the files for the output (except avp)
       
       
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
       open(unit=120, file='outflow.txt',
     & iostat=readerror)
                if (readerror.ne.0) then
                write(*,*) 'can''t open moments file'
                write(*,*) 'check path'
                stop
                end if
      open(unit=122, file='exchange.txt',
     & iostat=readerror)
                 if (readerror.ne.0) then
                 write(*,*) 'can''t open moments file'
                 write(*,*) 'check path'
                 stop
                 end if


       else fileopen

       ! generate the files for the output (except avp)
        
        
        open(unit=113, file='testcolenergy.txt',
     & iostat=readerror,position='append')
        if (readerror.ne.0) then
        write(*,*) 'can''t open moments file'
        write(*,*) 'check path'
        stop
        end if
        
         open(unit=101, file='disttrav.txt',
     & iostat=readerror,position='append')
        if (readerror.ne.0) then
          write(*,*) 'can''t open disttrav file'
          write(*,*) 'check path'
          stop
        end if
        open(unit=102, file='colcount.txt',
     & iostat=readerror,position='append')
        if (readerror.ne.0) then
          write(*,*) 'can''t open colcount file'
          write(*,*) 'check path'
          stop
        end if
        open(unit=103, file='timecounter.txt',
     & iostat=readerror,position='append')
        if (readerror.ne.0) then
          write(*,*) 'can''t open timecounter file'
          write(*,*) 'check path'
          stop
        end if
        open(unit=104, file='moments.txt',
     & iostat=readerror,position='append')
        if (readerror.ne.0) then
          write(*,*) 'can''t open moments file'
          write(*,*) 'check path'
          stop
        end if
      open(unit=120, file='outflow.txt',
     & iostat=readerror,position='append')
              if (readerror.ne.0) then
              write(*,*) 'can''t open moments file'
              write(*,*) 'check path'
              stop
              end if
      open(unit=122, file='exchange.txt',
     & iostat=readerror,position='append')
               if (readerror.ne.0) then
               write(*,*) 'can''t open moments file'
               write(*,*) 'check path'
               stop
               end if

      open(unit=121, file='momtot.txt',
     & iostat=readerror,action='read')
               if (readerror.ne.0) then
               write(*,*) 'can''t open moments file'
               write(*,*) 'check path'
               stop
               end if
      read(121,*) n_countav, mom_fiav, mom_secav, runcounter
      close(unit=121)
      end if fileopen



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
      testenergy=0.d0
      ! this is the loop over the particle type
      particletype: do pt=1, nop
      
      ! every particle can be run with different accuracy, this is the loop for that
      nparticleloop: do nparticle=1,np(pt)
     
     

      
      ! generate a particle 
      call generateparticle(pt)
      
      
      time=0.d0
      do while (parout==0)
     
      ! propagate the particle for a timeinterval
      call timestepprop(time)
      end do

      distsq=(particle(2)-disty/2.d0)**2+(particle(3)-distz/2.d0)**2
      
      if (distsq.le.radoutflow**2) then
        write(120,*) run, pt, nparticle, particle(1:6)
      end if
      
      

      
      end do nparticleloop
      runcounter(pt)=runcounter(pt)+dble(np(pt))
      end do particletype
      write(113,305) dble(run), testenergy
      CALL CPU_TIME(tnow)
      write(*,*) run, 'tnow',tnow
      if (floor(tnow/avpout)>floor(tprev/avpout)) then
        call printdens
        call printavp
        tprev=tnow
      end if
      
      if (mod(run,modmomout)==0) then
      write(101,300) disttrav
      volsubcell=product(distall)/dble(product(ncellall))
      do ip=1,nop
        colcount(ip,:,:,:,:)=colcount(ip,:,:,:,:)*
     &  influx(ip)/
     &  (dble(np(ip)*modmomout))
      colcount(ip,ip,:,:,:)=colcount(ip,ip,:,:,:)/2.d0
      end do
      write(102,300) colcount
      write(103,300) timecounter
      disttrav=0.d0
      colcount=0.d0
      timecounter=0.d0
      n_countav=n_countav+n_count3
      mom_fiav=mom_fiav+mom_fi3
      mom_secav=mom_secav+mom_sec3
        open(unit=121, file='momtot.txt',
     & iostat=readerror)
          if (readerror.ne.0) then
          write(*,*) 'can''t open moments file'
          write(*,*) 'check path'
          stop
         end if
        write(121,*) n_countav, mom_fiav, mom_secav, runcounter
        close(unit=121)
      write(104,*) n_count3, mom_fi3, mom_sec3
      n_count3=0.d0
      mom_fi3=0.d0
      mom_sec3=0.d0
    
      call flush
      end if
      if (mod(run,modnewavp)==0) then
      call newdens
      call newavp
      call flush
      end if
      
      
      end do
      call printdens
      call printavp
      n_countav=n_countav+n_count3
       mom_fiav=mom_fiav+mom_fi3
       mom_secav=mom_secav+mom_sec3
         open(unit=121, file='momtot.txt',
     & iostat=readerror)
           if (readerror.ne.0) then
           write(*,*) 'can''t open moments file'
           write(*,*) 'check path'
           stop
          end if
         write(121,*) n_countav, mom_fiav, mom_secav, runcounter
         close(unit=121)

      tprev=tnow
      ! close units and leave
      close(unit=101)
      close(unit=102)
      close(unit=103)
      close(unit=104)

      write(*,*) 'bye bye'
      
           
100   format(6A20)
200   format(10000000000I20)
300   format(10000000000E23.13)
305   format(10000000000E23.13)
306   format(10000000000E23.13)
      
      end program dsmc
