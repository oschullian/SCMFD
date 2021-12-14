      !
      !  Purpose: 
      !    This subroutine performs collisions between a particle and a background field
      !    
      !    
      !
      !
      !  Date       Programmer        Description of change
      !  ====       ==========        =====================
      !  04/05/14   Otto Schullian    
      !


      subroutine col
      use param
      implicit none
      
      ! the particle type of the field 
      integer :: pt2
      

      ! vectors needed to create random numbers (as soon as mkl library is available 
      ! this should be changed
      real, dimension(1) :: rvec
      real, dimension(3) :: rvece2
     
      ! type of the particle that collides with the field
      integer :: ptype
      
      ! loop integers
      integer :: k, l, m, n, s
      ! all initial/final absolute relative com velocities
      double precision, dimension(3) :: v1, v2, vcm, vr, vrf, v1f, v2f
      ! masses of the involved particle and the field
      double precision :: m1, m2, Mtot
      
      ! collision parameters: impact parameter, theta, phi
      double precision :: impar, theta, phi
      ! rotation matrices for the rotation of the final velocities and a orthogonal 
      ! transformation matrix that transforms vrel onto the x-axis
      double precision, dimension(3,3) :: Rot, rot2, mtmp
      ! three ortogonal vectors of mtmp
      double precision, dimension(3) :: e1, e2, e3
      
      ! index of the field-particle for the collision
      integer :: ind
      
      ! kinetic energies (for check)
      double precision :: en1, en2
      
      double precision, dimension(nop) :: kvalcdf
      integer :: ix, iy, iz
      
      integer :: check 
      
     
      ix=nint(particle(8))
      iy=nint(particle(9))
      iz=nint(particle(10))
      
      
      kvalcdf=kvalprev*ncelldens3(:,ix,iy,iz)
      do k=2,nop
      kvalcdf(k)=kvalcdf(k)+kvalcdf(k-1)
      end do
      kvalcdf=kvalcdf/kvalcdf(nop)
      
      
      
      call RANMAR(RVEC, 1)
      

      
      do k=1,nop
        if (rvec(1).le.kvalcdf(k)) then 
          pt2=k
          exit
        end if 
      end do 
      
      
      
      
      
      
      ! retrieve the type of particle that collides with the field
      ptype=nint(particle(7))
      
      
      ! collision is accepted, update collision counter
      if (switchCounter==1) then
      colcount(ptype,pt2)=colcount(ptype,pt2)+1.d0;
      else 
      switchCounter=1
      end if
 
      ! retrieve velocity of the particle that collides with the field
      v1=particle(4:6)
      

      ! call subroutine that generates v2, the velocity of the field particle that
      ! collides with the particle, this retrieve ind, the index of v2 within avp
      call generatecolvel(pt2,v2,ind)
      
      
      ! retrieve mass of each particle and calculate total mass
      m1=mass(ptype)
      m2=mass(pt2)
      Mtot=m1+m2
      
      ! calculate center of mass velocity and relative velocity
      vcm=1/Mtot*(m1*v1+m2*v2)
      vr=v1-v2
      
      check=0
      
      do while (check==0)
      call RANMAR(RVEC, 1)
      vrf(1)=(dble(rvec(1))-0.5d0)*2.d0
      call RANMAR(RVEC, 1)
      vrf(2)=(dble(rvec(1))-0.5d0)*2.d0
      call RANMAR(RVEC, 1)
      vrf(3)=(dble(rvec(1))-0.5d0)*2.d0
      
      check=1
      
      if (sum(vrf**2).gt.1.d0) check=0
      if (sum(vrf**2).le.0.d0) check=0
      
      
      
      end do 
      
    
      vrf=vrf*dsqrt(sum(vr**2))/dsqrt(sum(vrf**2))
      
      ! calculate the final velocities
      v1f=vcm+m2/Mtot*vrf
      v2f=vcm-m1/Mtot*vrf;
      
      ! calculate the initial and final kinetic energy
      en1=m1*sum(v1**2)+m2*sum(v2**2);
      en2=m1*sum(v1f**2)+m2*sum(v2f**2);
      
      ! compare if the energy is conserved
      testenergy(4)=testenergy(4)+en2-en1

          
      ! check for each component that the momentum is conserved
      do k=1,3
      testenergy(k)=testenergy(k)+m1*v1f(k)+m2*v2f(k)-m1*v1(k)-m2*v2(k)
      end do
      
      
      
      ! adjust avpvel, i.e. save the final velocity of the field particle 
      ! on the field vector
      particle(4:6)=v1f
      if (pt2==ptype) avpvel3(ind,pt2,:,ix,iy,iz)=v2f


      ! because v1f has changed, kval must be calculated always, therefore set
      ! kvalcalc to zero
      ! kvalcalc=0

      
      end subroutine col
      
      
      
      
