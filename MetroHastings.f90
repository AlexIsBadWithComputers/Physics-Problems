program metrohastings
  use random_phys581
  use statistics
  implicit none

  !use statistics
  integer::max,i,numbins
  integer,allocatable::binnies(:)
  real(kind = 8)::domain,std
  real(kind = 8),allocatable::points(:)

  max = 500000
  domain = 5.
  std =1.
  numbins = 150
  allocate(points(max))
  allocate(binnies(numbins))

  call MetroH(P,max,std,domain,points)

  call binValues(binnies,points,-1.0d0,1.0d0)
  open(12,file="points.txt",action='write')

  do i = 1,numbins
       write(12,*) binnies(i)
  end do



contains
  !this is the distribution fucntion
  real function P(x)
    implicit none
    real(kind=8), intent(in):: x

  !  P = 1./(2. * sqrt(3.14159)) * (sin(5.*x) + sin(2.*x) + 2.) * exp(-x**2)
    P = x**4 - x**2 + 0.1 * x
  end function P
  
  !this gives a unit gaussian distribution different standard deviation and
  !center.
  real function GaussChanger(Std,mu)
    implicit none
    real(kind=8),intent(in):: mu,std

    GaussChanger = random_normal() * std + mu

  end function GaussChanger

  !This first do-hickey shall compute the metropolis hastings algorithm for a given 
  !density P(x). I'll add more to this later when I understand it better. max is
  !the maximum interations and func is the density function to be passed in that bitch

  subroutine MetroH(func,max,std,domainsize,points)
    implicit none
    integer::i
    integer,intent(in):: max
    real(kind = 8), intent(in):: std,domainsize
    real(kind = 8):: y, xi,u,alpha
    real(kind = 8),intent(out),dimension(max)::points(:)
    real, external:: func

    call random_seed()
    !this is any number inside the domain
    call random_number(xi) 


    do i = 1, max
       !draw y = q(y|x)
     !  y = GaussChanger(std,xi)
       call random_number(y)
       y = 2.0 * y -1.
      ! print*, y
       alpha = min(1.0d0,func(y)/func(xi))


       call random_number(u)
       if (u.lt.alpha) then
          xi = y
          !continue algorithm
       else 
          xi = xi
          !stay here
       end if
       !save points
       points(i) = xi
    end do

  end subroutine MetroH



end program metrohastings
