      !
      !  Purpose: 
      !    This module shares data/parameters between subroutines.
      !    
      !  
      !  Code Summary: 
      !    The module defines parameters and arrays (with save statement)
      !    in order to exchange these between subroutines within the program.
      !
      !
      !  Date       Programmer        Description of change
      !  ====       ==========        =====================
      !  04/05/14   Otto Schullian    
      !

      
      module param

      
      !!!!!!!!!!!!!!
      ! constants  !
      !!!!!!!!!!!!!!
      double precision, parameter :: AK=0.831447147e4     !Boltzmann constant                                               !(in m**2 amu s**-2 K**-1)  
      double precision :: pi 
      double precision :: stpress=101325.e0
      double precision :: sttemp=273.15e0
      double precision :: AK2=1.3806488e-23 
      double precision :: mconv=1.66053892e-27            !conversion factor for amu in kg
      double precision :: clight=299792458.e0
      double precision :: planck=6.62606957e-34
      
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! random number generator parameters !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !seed for the random number generator
      integer :: seed
      ! parameters for the random number generator
      integer :: n1random=0
      integer :: n2random=0
      
      ! length of cumulative distribution function of a gaussian (i.e 
      ! the error function) - for the generation of the velocity on a wall
      integer, parameter :: ncdf = 100000
      ! velocity array of the cdf 
      double precision, dimension(ncdf) :: vcdf
      ! cdf of the errorfunction distribution function
      double precision, dimension(ncdf) :: cdf
      ! distance between two points on cdf
      double precision :: dcdf
      
      
      !!!!!!!!!!!!!!!!!!!!!!!!
      ! particle parameters  !
      !!!!!!!!!!!!!!!!!!!!!!!!
      ! number of particle types
      integer :: nop
      ! number of subtrajectories per full trajectory (i.e. if one particle 
      ! type needs to be simulated more often than another)
      integer, dimension(:), allocatable :: np
      ! reference particle diameter
      double precision, dimension(:), allocatable :: refdia
      ! reference temperature 
      double precision, dimension(:), allocatable :: reftemp
      ! viscosity index 
      double precision, dimension(:), allocatable :: visind
      ! mass 
      double precision, dimension(:), allocatable :: mass
      ! partial pressure (in Pa)
      double precision, dimension(:), allocatable :: pressure
      ! cross section pre factor (calculated with gammacalc)
      double precision, dimension(:,:), allocatable :: ccsfa
      ! average particle density 
      double precision, dimension(:), allocatable :: npartdens
      
      
      
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!
      ! background parameters  !
      !!!!!!!!!!!!!!!!!!!!!!!!!!
      ! number of particles in the background arrays
      integer :: noavp
      
      ! density of each particle species in every box and per particle type (first
      ! entry is particle type)
      double precision, dimension(:,:,:,:), allocatable :: ncelldens3
      ! vector that records the ncelldensity for a given trajectory, is used to update
      ! ncelldens3
      double precision, dimension(:,:,:,:),
     &  allocatable :: ncelldensnew3
      
      ! array containing the velocity distribution(s) of the background
      double precision, dimension(:,:,:,:,:,:), allocatable :: avpvel3
      ! array that saves the velocities for a trajectories and is later used to update 
      ! avpvel3
      double precision, dimension(:,:,:,:,:,:),
     & allocatable :: avpvelinter3
      
      
      ! counter to find the next position in array avpvelinter to save
      ! velocity
      integer, dimension(:,:,:,:), allocatable :: avtpos3
      
      ! counter that tells how many velocities were stored in avpvelinter3
      integer, dimension(:,:,:,:), allocatable :: avtexchange3
      
      ! counter of how many velocities are actually stored in avpvel3
      integer, dimension(:,:,:,:), allocatable :: nopcopos3
      
      ! counter that tells how many velocities in avpvel3 are maximally 
      ! exchanged in a updating procedure, i.e. the length of avpvelinter
      integer :: nvelchange
      
      
      
      !!!!!!!!!!!!!!!!!!!!!!
      ! system parameters  !
      !!!!!!!!!!!!!!!!!!!!!!
      
      ! sizes of the rectangular box in all directions
      double precision :: distx, disty, distz
      ! vector that contains all of the distances distx,disty,distz
      double precision, dimension(3) :: distall
      ! temperature of the wall(s)
      double precision :: tempwall
      ! number of cells in each dimension
      integer :: ncellx, ncelly, ncellz
      ! ncellx,ncelly,ncellz in one array
      integer, dimension(3) :: ncellall
      ! the change in pressure in z direction
      double precision :: dpdz
      ! the acceleration resulting from the pressure
      double precision :: accel
      ! maximal length of the timestep
      double precision :: dtmax
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!
      ! trajectory parameters  !
      !!!!!!!!!!!!!!!!!!!!!!!!!!
      
      ! length of a trajectory
      double precision :: trajlength
      ! number of trajectories run by the system
      integer :: totruns
      ! kval with every background gas type
      double precision, dimension(:), allocatable :: kvalprev
      ! particle vector, contains all information on particles
      double precision, dimension(10) :: particle
      ! integer that decides which avp and dens to read in
      integer :: avpstart
      ! if mod(no. trajectory,avpstart)==0 then the background is adjusted and
      ! avp is printed out
      integer :: modavpout
      ! counter for the avp that is printed out (used for names of files)
      integer :: nprint
      ! counter for the dens that is printed out (used for names of files)
      integer :: nprintdens
      ! percentage of velocities in avp that is either added, if avp has 
      ! not the full length yet, or exchanged if it is already full
      double precision :: inpcellperc
      

      !!!!!!!!!!!!!!!!!!!!!
      ! output parameters !
      !!!!!!!!!!!!!!!!!!!!!
      
      ! random offset for dtprint (for every particle), so that there is not 
      ! bias (for example if they all were registered at t=0)
      double precision :: dtprintoffset
      ! the timeintervals for saving the velocity into avpvelinter 
      double precision, dimension(:), allocatable :: dtprint
      ! a switch that turns on the registration of the colcount, disttrav,
      ! and timecounter
      integer :: switchCounter
      ! time for printing avp
      double precision :: tprev, tnow
      
      
      !!!!!!!!!!!!!!!
      ! output data !
      !!!!!!!!!!!!!!!
      
      ! colcount(p1,p2) number of collisions of particle p1 with 
      ! background p2
      double precision, dimension(:,:), allocatable :: colcount
      ! distance travelled of a particle (type)
      double precision, dimension(:), allocatable :: disttrav
      ! total time the particle (type) has traveled
      double precision, dimension(:), allocatable :: timecounter
      
      ! contains the cumulative error of the pre and postcollisional
      ! momenta and energy
      double precision, dimension(4) :: testenergy
      
      ! counts the number of times a particle has hit the wall, first index
      ! is which wall (x,y,z) and second is whether it hit the left (*,1) or
      ! right (*,2) side
      integer, dimension(3,2) :: wallhitstat
      
      ! arrays for the moments
      ! normalisation for the moments
      double precision, dimension(:,:,:,:), allocatable :: n_count3
      ! sum for first moment, i.e. first moment=mom_fi3/n_count
      double precision, dimension(:,:,:,:,:), allocatable :: mom_fi3
      ! sum for second moment, i.e. second moment=mom_sec3/n_count
      double precision, dimension(:,:,:,:,:), allocatable :: mom_sec3
      ! normalisation for the moments
      double precision, dimension(:,:,:,:), allocatable :: n_countav
      ! sum for first moment, i.e. first moment=mom_fi3/n_count
      double precision, dimension(:,:,:,:,:), allocatable :: mom_fiav
      ! sum for second moment, i.e. second moment=mom_sec3/n_count
      double precision, dimension(:,:,:,:,:), allocatable :: mom_secav

      

     
           
      end module param
