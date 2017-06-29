!This program contains a function which will calculate the energy of photons
!coming through a slab after they have been binned appropriately.
program intensity
  use statistics
  implicit none 
  double precision:: zmax, prob,min,max
  integer, allocatable:: taubin(:),mubin(:),phibin(:)
  double precision,allocatable:: tauarr(:),muarr(:),phiarr(:),inten(:),vals(:),junk(:)
  integer::binlength,i,rvarnum,NPhotons,binsize,j,k
  NPhotons = 10**6
  zmax = 10.
  prob = 1.
  binsize = 20
  min = 0.
  max = 90.
  allocate(inten(binsize))
  allocate(mubin(binsize))
  allocate(vals(binsize))
  allocate(junk(NPhotons))
  !set value of each bin
  do i = 1,binsize
     vals(i) = acos((1.*i)/20.)!)*180./3.1415962
  end do
  !print*, vals
  call PhotonScatter(NPhotons,muarr,phiarr,zmax,prob)
  muarr = acos(muarr)*180./3.1415962
  call binValues(mubin, muarr,min,max)

 !MAKE SURE YOU YOU REVERSE THE MUBIN ARRAY BECAUSE IT IS BACKWARDS ONCE YOU GET ANGLES.
  !call intensityCalc(mubin,vals,Inten,20,NPhotons)
  !print*, minval(muarr),maxval(muarr),mubin
  open(12,file="radiativetransfer.txt",action="write")
  
 ! print*,minval(junk),maxval(junk), "distingishe me more"

close(12)

contains 
  !This subroutine calculates the normalized intensity of data
  !of photons scattering through a slab. All arrays must be
  !double precision to prevent computer saddness. 
  subroutine IntensityCalc(NumBin,Angle,Intensity,Size,NPhotons)
    implicit none
    integer, intent(in):: Size, NPhotons
    integer,dimension(size),intent(in):: NumBin
    double precision, dimension(Size), intent(in):: Angle
    double precision, dimension(Size), intent(out):: Intensity
    integer:: i

    do i = 1, Size
       Intensity(i) = (NumBin(i)*1.0* Size)/(2. * NPhotons * Angle(i))
    end do

    return
  end subroutine IntensityCalc
  !This makes a random array for the tau values modified to return only one number
  subroutine tau(tauu)
    implicit none
    !integer,intent(in):: size
   ! doubleprecision,intent(out),dimension(size)::tauarray
    doubleprecision,intent(out)::tauu
    doubleprecision::harvest
    integer::i

    call random_seed()
    !do i =1,size
       call random_number(harvest)
       tauu = -log(1-harvest)
    !end do
  end subroutine tau

  !random array of theta bvalues
  subroutine mu(muu)
    implicit none
    !integer,intent(in):: size
    doubleprecision,intent(out)::muu
    doubleprecision::harvest
    integer::i

    call random_seed()

    !do i =1,size
       call random_number(harvest)
       muu  = 2*harvest - 1
   ! end do

  end subroutine mu
  !random array of phi values
  subroutine phi(phia)
    implicit none
   ! integer,intent(in):: size
    doubleprecision,intent(out)::phia
    doubleprecision::harvest,pi
    integer::i
    pi = 3.14159
    call random_seed()

  !  do i =1,size
       call random_number(harvest)
       phia =2*pi* harvest
  !  end do


  end subroutine phi

  !This one simulates scattering of a photon through a slab
  subroutine PhotonScatter(NPhotons,packmu,packphi,zmax,prob)
    implicit none
    integer, intent(in)::NPhotons
    double precision,intent(in):: zmax,prob
    double precision, dimension(NPhotons):: muarray, phiarray
    double precision,allocatable,intent(out):: packmu(:),packphi(:)
    logical, dimension(NPhotons):: fudge
    integer::i
    double precision:: x,y,z,AbsOrRef,tauu,muu,phii,Rmu,Rtau,Rphi,pi,a_or_s
    !it is important to assume that our slab is infinite in the x
    !and y direction for simplicity.
    pi = 3.14159

    !initial positions - always origin

  
    !To begin, we need to sample phi from 0 - 2pi and theta from the
    !function mu  =cos(theta)= 2*random - 1 to calculate mu.
    !as well tau = -ln(1-random)
    do i = 1, NPhotons
       x = 0.
       y = 0.
       z = 0.
       inslab: do while (z.ge.0 .and. z.lt.zmax)
         call mu(muu)
         call phi(phii)
         call tau(tauu) 

        !  tau = tauarray(i)
        !  mu = 2.*Rmu - 1 
        !  phi = 2.*pi*Rphi


          x = x + tauu*sin(phii)*muu
          y = y + tauu*sin(acos(muu))*sin(phii)
          z = z + tauu*muu

          if (z.lt.0) then
             !photon is reflected :(
             fudge(i) = .false.
          end if
          call random_number(a_or_s)
          !if reflective only, this loop should never end
          if (prob .lt. a_or_s) then
             fudge(i) = .false.
             print*, "ramalamadingdong"
             !photon is absorbed. No longer does anything.
             exit
          else if(z .gt. zmax) then
             muarray(i) = muu
             phiarray(i) = phii
             fudge(i) = .true.
             exit
          end if
       end do inslab
    end do
  
    print*,maxval(muarray), minval(muarray), "distinguish me"
    packmu = pack(muarray,fudge)
    packphi = pack(phiarray,fudge)
    print*, size(fudge) - size(packmu)

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
