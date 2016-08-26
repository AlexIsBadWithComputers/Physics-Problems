program sa
  use statistics
  implicit none
  integer:: max,i,junk,j,max2
  real(kind=8), dimension(:),allocatable::points
  integer,dimension(:),allocatable:: bins
  
  max = 1000
  max2 = 100
  junk = 0
  allocate(bins(100))

 ! allocate(points(max))
  call simulatedannealing(fv,fv,points,max,max2)
  call binvalues(bins,points,-1.d0,1.d0) 
  
  open(15,file="siman3.txt")
  do i = 1, size(bins)
     write(15,*) bins(i), (2.0*i/100 - 1.) - 1./50
  end do
  write(15,*) ""
  write(15,*) ""
  do i = 1,size(points)
     write(15,*)i, points(i)
  end do
  print*, average(points,size(points)), "look at me"
contains
!taks an external function to find minima of func, a cost function cost
!as well requires an allocatable rank 1 array points to return
!your values. max1 and max2 are the number of interations at a temperature
!and the number of temperature decreaces respectivbely. Currently, the
!temperature schedule is very quick, as it simply scales as 1/j, hwere
!j = 1,max2


  real function fv(x)
    implicit none
    real, intent(in):: x
    fv = x**4 - x**2 + 0.1 * x
  !  print*, "smd", fv
  end function fv
  
  real function tschedule(x)
    implicit none
    real(kind = 8), intent(in):: x
    tschedule = .1
  end function tschedule
end program sa
