!This program contains a function which will calculate the energy of photons
!coming through a slab after they have been binned appropriately.
program intensity
  use statistics
  implicit none 
  double precision:: zmax, prob,min,max,NPhotons
  integer, allocatable:: taubin(:),mubin(:),phibin(:),switch(:)
  double precision,allocatable:: tauarr(:),muarr(:),phiarr(:),inten(:),vals(:),bvals(:)
  integer::binlength,i,rvarnum,binsize,j,k
  NPhotons = 10**6
  zmax = 10.
  prob = .5 !scattering probability
  binsize = 20
  min = 0.
  max = 90.
  allocate(inten(binsize))
  allocate(mubin(binsize))
  allocate(vals(binsize))
  allocate(switch(binsize))
  allocate(bvals(binsize))
  
  do i = 1,binsize
     vals(i) = 0.99999999999*i/20. - .99999999999/40.
  end do
  
  call PhotonScatter(NPhotons,muarr,phiarr,zmax,prob)
  
!print*, muarr
!  
 
  call binValues(mubin,muarr)
  
  
  print*, mubin
  print*,
  print*,
    

 ! print*, mubin
  call intensityCalc(mubin,vals,Inten,20,1.0*sum(mubin))
!
  !print*, minval(muarr),maxval(muarr),mubin
!  open(12,file="radiativetransferabsorby.txt",action="write")
 
!  do i = 1, binsize
     
 !    write(12,*) inten(i), acos(vals(i))*180./3.1415962, mubin(i)
     
 ! end do

 !close(12)
contains 
  !This subroutine calculates the normalized intensity of data
  !of photons scattering through a slab. All arrays must be
  !double precision to prevent computer saddness. 
  subroutine IntensityCalc(NumBin,muval,Intensity,Size,NPhotons)
    implicit none
    integer, intent(in):: Size
    real,intent(in)::NPhotons
    integer,dimension(size),intent(in):: NumBin
    double precision, dimension(Size), intent(in):: muval
    double precision, dimension(Size), intent(out):: Intensity
    double precision:: num, dom
    integer:: i

    do i = 1, Size
       num = NumBin(i) * 1.0 * Size
       dom = 2. * NPhotons * muval(i)! (acos(muval(i))*180/3.14159)
       intensity(i) = num / dom 
       print*, intensity(i), num,dom, numbin(i)

    end do 
    return
  end subroutine IntensityCalc
  !This makes a random array for the tau values
  subroutine tau(size,tauarray)
    implicit none
    integer,intent(in):: size
    doubleprecision,intent(out),dimension(size)::tauarray
    doubleprecision::harvest
    integer::i

    call random_seed()
    do i =1,size
       call random_number(harvest)
       tauarray(i) = -log(1-harvest)
    end do
  end subroutine tau

  !random array of theta bvalues
  subroutine mu(size,thetaarray)
    implicit none
    integer,intent(in):: size
    doubleprecision,intent(out),dimension(size)::thetaarray
    doubleprecision::harvest
    integer::i

    call random_seed()

    do i =1,size
       call random_number(harvest)
       thetaarray(i) = -acos(harvest )
    end do

  end subroutine mu
  !random array of phi values
  subroutine phi(size,phiarray)
    implicit none
    integer,intent(in):: size
    doubleprecision,intent(out),dimension(size)::phiarray
    doubleprecision::harvest,pi
    integer::i
    pi = 3.14159
    call random_seed()

    do i =1,size
       call random_number(harvest)
       phiarray(i) = harvest
    end do


  end subroutine phi

  !This one simulates scattering of a photon through a slab
  subroutine PhotonScatter(NPhotons,packmu,packphi,zmax,prob)
    implicit none
    double precision,intent(in):: zmax,prob,NPhotons
    double precision, dimension(int(NPhotons)):: muarray, phiarray
    double precision,allocatable,intent(out):: packmu(:),packphi(:)
    double precision,allocatable::trackypoo(:,:)
    logical, dimension(int(NPhotons)):: fudge
    integer::i,j,count,kill
    double precision:: x,y,z,AbsOrRef,tau,mu,phi,Rmu,Rtau,Rphi,pi,a_or_s,theta
    !it is important to assume that our slab is infinite in the x
    !and y direction for simplicity.
    pi = 3.14159
    allocate(trackypoo(100000,3))

    !initial positions - always origin
    kill = 0
    call random_seed()
    !To begin, we need to sample phi from 0 - 2pi and theta from the
    !function mu  =cos(theta)= 2*random - 1 to calculate mu.
    !as well tau = -ln(1-random)
    do i = 1,int(NPhotons)
       x = 0.
       y = 0.
       z = 0.
       count = 0
       inslab: do while (z.ge.0 .and. z.lt.zmax)
          call random_number(Rmu)
          call random_number(Rphi)
          call random_number(Rtau)
          call random_number(a_or_s)
          
          
          tau = -log(Rtau)!-log(1.- Rtau)
          mu = 2.*Rmu - 1. 
         ! theta = -acos(2.*Rmu-1.)
          phi = 2.*pi*Rphi
          
          count = count + 1
         
          
          x = x + tau/zmax * sin(acos(mu)) * cos(phi)!mu
          y = y + tau/zmax * sin(acos(mu)) * sin(phi)
          z = z + tau/zmax * mu
          
          trackypoo(count,1) = x
          trackypoo(count,2) = y
          trackypoo(count,3) = z

          if (z.lt.0) then
             !photon is reflected :(
            ! print*, ":-("
             muarray(i) = cos(theta)!mu
             fudge(i) = .false.
             exit
          end if
          if(z .gt. zmax) then
           !  print*,"Hello!"
             muarray(i) =  mu
             phiarray(i) = phi
             fudge(i) = .true.
             open(34,file = "path.txt",action="write")
             
             do j = 1,count
                write(34,*) trackypoo(j,:)
             end do
             kill = 1
             close(34)
             exit
          end if
                   
     ! if (prob .gt. a_or_s) then
        !   print*,prob, a_or_s
      !     fudge(i) = .false.
           !photon is absorbed. No longer does anything.
       !    exit
        ! end if
          
       end do inslab
       if (kill == 1) exit
       trackypoo = 0
    end do
  !  print*,maxval(muarray), minval(muarray), "distinguish me"
    packmu = pack(muarray,fudge)
    packphi = pack(phiarray,fudge)
    print*, size(packmu)
    deallocate(trackypoo)
  end subroutine PhotonScatter


end program intensity
!  rvarnum = 10**7   
  !call tau(rvarnum,tauarr)
!call mu(rvarnum,muarr)
!call phi(rvarnum,phiarr)

!do i = 1,rvarnum
!   muarr(i) = acos(muarr(i))
!end do
!binlength = 100
!allocate(taubin(binlength))
!allocate(mubin(binlength))
!allocate(phibin(binlength))
!print*, "hey"


!call binvalues(taubin,binlength,tauarr)

!call binvalues(mubin,binlength,muarr)
!call binvalues(phibin,binlength,phiarr)

!open(12,file="histogram.txt",action="write")

!do i =1,binlength
!   write(12,*) taubin(i), i

!end do
!write(12,*) ""
!write(12,*) ""
!print*, "one"
!do i =1,binlength
!   write(12,*) mubin(i), i
!end do

!write(12,*) ""
!write(12,*) ""
!print*, "two"
!do i =1,binlength
!   write(12,*) phibin(i), i
!end do

!write(12,*) ""
!write(12,*) ""


!close(12)
