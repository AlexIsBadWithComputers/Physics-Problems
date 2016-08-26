module statistics 
  !Module containing subroutines and functions for statistics related
  !calculations. Alex + Zach contributuinos. 
  use random_phys581
  implicit none 

contains 

  !Bins values over the closed interval [LOWER:UPPER] 
  !which are optional arguments which default to [0:1] 
  !The number of bins is defined by the size of the array
  !provided as an arggument. The binned vales are input 
  !as an array of size numvalues and the number of entries 
  !per bin is returned as an integer vector.
  !Values outside the lower and upper bounds are discarded 
  !and a warning is printed. 
  subroutine binValues(BINS, VALUES, LOWER, UPPER)
    implicit none
    double precision, OPTIONAL :: LOWER, UPPER 
    integer, DIMENSION(:), INTENT(OUT) :: BINS
    double precision, DIMENSION(:), INTENT(IN) :: VALUES
    integer:: i, j, num_bins
    double precision:: x
    double precision:: low=0., high=1., range=1.0

    num_bins = size(BINS)

    if (present(LOWER) ) then 
       low = lower
    end if
    if (present(UPPER) ) then
       high = upper
    end if
    
    range = high - low
    
    BINS = 0
    do j=1, size(VALUES)
       x = VALUES(j)
       bin: do i=1,num_bins
          !Warn if values given are not within the provided bounds. These values are discarded.
          if (x .lt. low ) then
             write(*,*) "In binning sub, Value x: ", x, " is lower than lower binning limit ", low
             exit bin
          end if
          if (x .gt. high) then
             write(*,*) "In binning sub, Value x: ", x, " is higher than upper binning limit", high
             exit bin
          end if
             
          if ( x .le. low + range*i/real(num_bins) ) then
             BINS(i) = BINS(i)+1
             exit bin
          end if
       end do bin
    end do

  end subroutine binValues


  !Calculates the table chi square value given a set 
  !of observed and expected values in arrays of equal size.
  subroutine chiSquaredTest(observed, expected, chiSquared) 
    implicit none
    double precision, INTENT(IN):: observed(:), expected(:)
    double precision, INTENT(OUT):: chiSquared
    integer :: i 

    if (size(observed).ne. size(expected) ) then 
       print*, "In Chi square test, observed and expected arrays size mismatch"
    end if 
    
    chiSquared = 0.0
    do i=1, size(observed)
       !write(*,*) chiSquared, observed(i), (observed(i) - expected(i))**2/expected(i)
       chiSquared = chiSquared + (observed(i) - expected(i))**2 / expected(i)
    end do

  end subroutine chiSquaredTest

  !Given an existing mean of N numbers, the new mean 
  !of N+1 values where x is the new value to add to the 
  !mean is returned by this function.
  pure double precision function addToMean(x, mean, N)
    double precision, INTENT(IN):: x, mean
    integer, INTENT(IN) ::N
    
    addToMean = mean + (x-mean)/(N+1)
  end function addToMean

  !Given an existing mean of N numbers, the new mean 
  !if we remove a value x from it is returned by this 
  !function.
  pure double precision function subFromMean(x, mean, N)
    double precision, INTENT(IN) :: x, mean
    integer, INTENT(IN)::N

    subFromMean = (1.0*N*mean-x)/(N-1)
  end function subFromMean

  !Given a function (pointer) p and an input array xvalues of size N, 
  !the subroutine will run a random-walk metropolis-hastings algorithm 
  !with a gaussian proposal distribution. 
  !The optional arguments x0 and sigam affect the starting position
  !and the width of the proposal gaussian. 
  !The optional arrays acceptRate and acceptProb will be returned with 
  !the accptance rate and acceptance probability at each step. 
  subroutine metropolis(p, xvalues, N, x0, sigma, acceptRate, acceptProb)
    implicit none 
    double precision, external :: p
    double precision, INTENT(OUT) :: xvalues(N), acceptRate(N), acceptProb(N)
    double precision, INTENT(IN), OPTIONAL :: x0, sigma
    integer, INTENT(IN) :: N
    OPTIONAL :: acceptRate, acceptProb
    integer :: i
    double precision :: x, sig, y, u, alpha, mean, avgAccept=0.0

    x = 0.0
    if (PRESENT(x0)) x = x0

    sig = 1.0
    if (PRESENT(sigma)) sig = sigma

    do i=1, N
       !Draw y ~ q(y|x)
       mean = x
       y = sig*random_normal() + mean
       !a(x, y) <- Min(1, p(y)q(y|x) / p(x)q(x|y) )
       alpha = min(1.0, p(y)/p(x)*1.0) !One for symmetric distribution
       !Draw u~Unif(0,1)
       call random_number(u)
       !If u < a(x,y) then Accept:
       if (u .lt. alpha) then
          avgAccept = addToMean(1.0d0, avgAccept, i)
          x = y !Accept as move
       else
          avgAccept = addToMean(0.0d0, avgAccept, i)
          x = x !Stay
       end if

       xvalues(i) = x
       if( PRESENT(acceptRate)) acceptRate(i) = avgAccept
       if( PRESENT(acceptProb)) acceptProb(i) = alpha

    end do

  end subroutine metropolis



!This function takes the average of a 1D array, length N
     real(kind = 8) function Average(array,N) 
      implicit none
      integer,intent(in):: N
      real(kind = 8),intent(in), dimension(N):: array
      
      Average = sum(array(1:N))/(1.0*N)
     ! print*, "I DONE DID GOOD"
    
    end function Average
    !This function takes a 1D array of length N of random numbers and calculate
    !s the autocorrelation function through a shift of k values.
    real(kind =8) function AutoCorr(RanNums,N,k)
      implicit none
      integer,intent(IN)::N,k
      real(kind = 8),intent(in), dimension(N)::RanNums
      real(kind = 8):: xbar,numerator, denominator
      integer:: t
    
      xbar = Average(RanNums,N) !find average of random numebers
     
      numerator = 0             !initialize sums
      denominator = 0 
     
      if (N-k .lt. 1) then
        print*, "You don't know what you're doing"
     else
        do t = 1, N - k           !apply sumation to equation (11) 
           numerator = numerator + (RanNums(t) - xbar)*(RanNums(t + k) - xbar)
           denominator = denominator + (RanNums(t) - xbar)**2
        end do
     
      AutoCorr = abs(numerator/denominator)
     end if 
    end function AutoCorr
    
    !calculates teh standard deviation of a discrete set returning the 
    !standard deviation sig, from the set called set .
    real(kind = 8) function StandDev(set)
      implicit none
      real(kind = 8),dimension(:):: set
      real(kind = 8):: mu,sum,ave
      integer::N,i
      
      N = size(set)
      ave = Average(set,N)
      
     set = set - ave
     sum = 0
     do i = 1,N
        sum = sum + set(i)**2
     end do
     
     StandDev = sqrt((1./ave) * sum)
   end function StandDev



    !this function runs the accept reject method on an external function p,
  !and takes input for number of interations n, maximum value pmax, 
  !as well as range on which p is constrained a and b. Data is returned
  !as an n dimensional array data.
   subroutine AcceptReject(p,pmax,a,b,n,data)
     implicit none
     real,external::p
     real(kind = 8),intent(in):: pmax,a,b
     integer,intent(in)::n
     real(kind=8),dimension(n),intent(out):: data
     logical:: accept
     real(kind = 8):: try,reject,xtry
     integer::i

     call random_seed()

     do i=1,n
        accept =.false. !only write values onc they are accepted.
        do while(accept .eqv. .false.)
           call random_number(try)
           call random_number(reject)

           xtry = a + (b - a) * try !xtry uniformly distributed over [a,b)

           if (p(xtry) .ge. reject*pmax) then !accept 
              data(i) = xtry
              accept = .true.
              exit
           end if
        end do
     end do
   end subroutine AcceptReject

   !same as above, but someone was whiney and wanted a function too.
   !remember to call random_seed() before use, or else shit wont work.
    real(kind = 8) function AcceptRejectFunc(p,pmax,a,b)
     implicit none
     real(kind = 8),external::p
     real(kind = 8),intent(in):: pmax,a,b
     logical:: accept
     real(kind = 8):: try,reject,xtry
     integer::i

        accept =.false. !only write values onc they are accepted.
        do while(accept .eqv. .false.)
           call random_number(try)
           call random_number(reject)

           xtry = a + (b - a) * try !xtry uniformly distributed over [a,b)

           if (p(xtry) .ge. reject*pmax) then !accept 
              AcceptRejectFunc = xtry
              accept = .true.
              exit
           end if
        end do

   end function AcceptRejectFunc



!taks an external function to find minima of func, a cost function cost
!as well requires an allocatable rank 1 array points to return
!your values. max1 and max2 are the number of interations at a temperature
!and the number of temperature decreaces respectivbely. Currently, the
!temperature schedule is very quick, as it simply scales as 1/j, hwere
!j = 1,max2
subroutine SimulatedAnnealing(func,cost,points,max1,max2)
    implicit none
   integer::i
    integer,intent(in):: max1,max2!
    !real,intent(in)::tstar!t
    real:: u,alpha,c_new,c_old,dt,T,tfin,cst
    real(kind = 8),intent(out),dimension(:),allocatable::points(:)
    real, external:: func,cost
    integer::count,j   
    allocate(points(max1*max2))
    call random_seed()
    !this is any number inside the domain,gets things rolling
    call random_number(c_old) 
    c_old = 2. * c_old - 1.
    dt = 1.
    count = 0
    t = 1.
    print*, c_old,func(c_old),t
  do j = 1,max2
    do i = 1, max1
      
       call random_number(c_new)
       c_new = 2. * c_new - 1.
       cst = cost(c_new) - cost(c_old)
   ! 
       if (cst.le.0.) then
          c_old = c_new
 
       else
          alpha = exp(-cst/T)
          call random_number(u)
          if (u.le.alpha) then
             c_old = c_new
       
          end if
       end if
      !save points
       count = count + 1
    !   print*, count,t
       points(count) = c_old
!       
!       
    end do
    !count = count + 1
    T= 1.0/real(j)
  end do
    
  end subroutine SimulatedAnnealing

 end module statistics

