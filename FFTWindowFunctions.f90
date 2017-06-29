module fftstuff
implicit none
!This is usefull stuff for FFT. Mostly windows. Look them up if you've forgotten. 
contains
subroutine hanningWindow(vec,window) 
  implicit none
  complex(kind = 8), intent(in),dimension(:) :: vec
  complex(kind =8), intent(out),dimension(size(vec)):: window
  integer::i
  
  
  do i = 1, size(vec) 
     window(i) =(0.5d0 - 0.5d0 * cos(2.d0*3.1415962d0 * real(i,kind =8)/(size(vec)-1)))*vec(i)
  
  end do
end subroutine hanningwindow

subroutine hammingwindow(vec,window)
  implicit none
  complex(kind = 8), intent(in),dimension(:)::vec
  complex(kind = 8), intent(out),dimension(size(vec))::window
  integer:: i

  do i = 1, size(vec)
     window(i) = (0.54d0 - 0.46d0*cos(2.d0*3.1415962d0*real(i,kind = 8)/(size(vec,kind=8) - 1.d0)))*vec(i)
  end do
end subroutine hammingwindow

subroutine blackmann(vec,window)
  implicit none
  complex(kind = 8), intent(in),dimension(:)::vec
  complex(kind = 8), intent(out),dimension(size(vec))::window
  real(kind = 8):: pi
  complex(kind = 8):: factor
  integer:: i
  
  pi = 3.1415962
  
  do i = 1, size(vec)
     factor = 0.42d0-0.5d0*cos(2.d0*pi*real(i,kind =8)/(size(vec,kind = 8)-1.d0)) &
               +0.08d0*cos(4.d0*pi*real(i,kind=8)/(size(vec,kind=8)-1.d0))
     window(i) = factor * vec(i)
  end do
end subroutine blackmann

subroutine BlackHarris(vec,window)
  implicit none
  complex(kind = 8), intent(in),dimension(:)::vec
  complex(kind = 8), intent(out),dimension(size(vec))::window
  real(kind = 8):: pi
  complex(kind = 8):: factor
  integer:: i
 
  pi = 3.1415962

  do i = 1, size(vec)
     factor = 0.35875d0 - 0.48829d0 *cos(2.d0*pi * real(i,kind=8)/(size(vec,kind=8) - 1.d0)) &
              + 0.14128d0 * cos(4.d0 * pi * real(i,kind=8)/(size(vec,kind=8) - 1.d0)) &
              - 0.01168d0 * cos(6.d0 * pi * real(i,kind=8)/(size(vec,kind=8) - 1.d0))
     
     window(i) = factor * vec(i)
  end do
end subroutine BlackHarris

subroutine fftflip(vec)
  implicit none
  complex(kind = 8), dimension(:),intent(inout)::vec
  complex(kind = 8), dimension(size(vec)/2)::flipy,flipy2
  integer:: N 
  N = size(vec)
  
  flipy = vec(1:N/2)
  flipy2 = vec(N/2+1:N)
  vec(N/2+1:N) = flipy
  vec(1:N/2) = flipy2
end subroutine fftflip

subroutine fftflip2(vec)
  implicit none
  real, dimension(:),intent(inout)::vec
  real, dimension(size(vec)/2)::flipy,flipy2
  integer:: N 
  N = size(vec)
  
  flipy = vec(1:N/2)
  flipy2 = vec(N/2+1:N)
  vec(N/2+1:N) = flipy
  vec(1:N/2) = flipy2
end subroutine fftflip2

real(kind=8) function timestep(N,dt,i)
  implicit none
  integer,intent(in)::N,i
  real(kind = 8)::dt

  timestep = -real(N,kind=8)/(2.d0 * (dt)) + real(i,kind = 8)/dt
end function



subroutine filelength(filename, length)
  implicit none
  character(len = *),intent(in):: filename
  integer,intent(out):: length
  character(len=1):: junk
  integer:: max, j, ios
  open(12, file = filename, action = 'read')
  print*, filename
  length = 0
  max = 100000000
 

  do j = 1,max
     read(12,*,IOSTAT=ios) junk
    if (ios .ne. 0) then
        exit
    end if
    if (j == max) then
        write(*,*) "Too many dicks on the dance floor"
        write(*,*) "(Too many dicks)"
        write(*,*) "But for real, file has more than 100000000lines"
        stop
    end if
    length = length + 1
  !  print*, length
 !  print*,lines
  end do
 ! print*, length
  close(12)  
end subroutine filelength


integer function poweroftwo(number)
  implicit none
  integer::n_temp,n
  integer,intent(in)::number

  n_temp = log(real(number))/log(2.)
  n = ceiling(real(n_temp)) !truncate data to nearest power of two
  poweroftwo = 2**n
end function poweroftwo
end module
