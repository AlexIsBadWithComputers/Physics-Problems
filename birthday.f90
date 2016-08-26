program birthday
implicit none
integer:: i,homies,trials,j,k,success
integer,allocatable:: dat(:)
real::draw
!classic birthday problem
success = 0
homies = 30
trials = 1000000
allocate(dat(homies))
call random_seed()
big: do k = 0, trials
       do i = 1,homies
          call random_number(draw)
          dat(i) = int(draw * 364)
       end do
  
   !Note: the first element of the array is compared
   !with each element besides itself, then as we make our way through 
   !the loop we no longer need to compare that element with any other elements
   !of the array, hence the i+1 dependancy on j.
       do i =1,size(dat)
          do j = i+1,size(dat) !ignore elements that have already been counted
             if (dat(i) == dat(j)) then
                success = success + 1
                cycle big !Only looking for single match
             end if
          end do
       end do
    end do big

print*, 1.0*success/(trials), success

 



end program birthday
