program ising
!This computes an ising model for a theoretical magnet slab
!at various temperatures. I didn't comment very well. sorry.
!but it uses metropolis-hastings to decide if the spins flip or not
!Although I like to call it one temperature simulated annealing because
!i'm a rule breaker
implicit none
integer::siz,sweeps,i,j
real(kind = 8),allocatable,dimension(:)::en,men,magn,mmagn
real(kind = 8),dimension(4)::temp

temp = (/3.0,2.5,2.25,1.0/)
siz = 32 
sweeps = 100000 
open(16, file = "isingtest.txt")

allocate(en(sweeps),men(sweeps),magn(sweeps),mmagn(sweeps))
!call simulation(temp(4),sweeps,siz,en,men,magn,mmagn)

do i = 1, 4
  print*, i
  en =0. ; men = 0. ; magn =0; mmagn= 0.;
  call simulation(temp(i),sweeps,siz,en,men,magn,mmagn)
  do j = 1,sweeps
    write(16,*) j, en(j),men(j),magn(j),mmagn(j)
   ! print*,  j, en(j),men(j),magn(j),mmagn(j)
  end do
  print*,mmagn(size(mmagn)), "MMAG @: ",temp(i)

  write(16,*)""
  write(16,*)""
end do
close(16)
deallocate(en,men,magn,mmagn)

contains
subroutine simulation(temp,sweeps,siz,en,meanen,magnet,meanmagnet)
  implicit none
  real(kind=8),intent(in)::temp
  integer,intent(in):: sweeps,siz
  real(kind = 8),dimension(siz,siz)::spin
  real(kind = 8),dimension(sweeps),intent(out)::en,meanen,magnet,meanmagnet
  integer:: i,j,k,l,flipy
  real(kind = 8)::r1,r2,r,d
  spin= +1.

  call random_seed()
  
  en(1) = ToToNRG(spin)
  meanen(1) = en(1)
  magnet(1) = magnetizm(spin)
  meanmagnet(1) = abs(magnet(1))
  flipy = 0
  
  do l = 2,sweeps
    do k = 1, siz**2 - 1
        call random_number(r1)
        call random_number(r2)

        i = int(r1 * siz ) + 1  !Need to generate numbers between [1,siz]
        j = int(r2 * siz ) + 1
        d = -4.0 * delt(spin,i,j)
        if (d.lt.0) then
          spin(i,j)=-1.0 * spin(i,j)
          flipy = flipy + 1
        else
        call random_number(r)
          if (exp(-d/temp ).gt. r) then 
            spin(i,j) = -1.0 * spin(i,j)
            flipy = flipy + 1
          end if
        end if
    end do
 
    en(l) = totonrg(spin)
    meanen(l) = (meanen(l-1) * (l - 1) + en(l))/real(l,kind = 8)
    magnet(l) = magnetizm(spin)
    meanmagnet(l) = (meanmagnet(l-1) * (l-1) + magnet(l))/real(l,kind = 8)
    !print*,l, en(l), meanen(l), magnet(l), meanmagnet(l)
    !write(16,*) l, en(l), meanen(l), magnet(l)
  end do
!close(16)
  print*, "Flips: ",flipy

end subroutine simulation

!Periodic boundaries shake and bake. numpnts + 1, 0 becomes numpts and so on
real(kind = 8) function boundaries(s,i,j)
  implicit none
  real(kind = 8),intent(in),dimension(:,:):: s
  integer, intent(in):: i,j
  integer::xind,yind
 
  !assign periodic spins all lazy like
  xind = modulo((i - 1), size(s,2)) + 1
  yind = modulo((j - 1), size(s,2)) + 1
  boundaries = s(xind,yind)
end function boundaries

real(kind = 8) function delt(s,i,j)
  implicit none
  real(kind = 8),dimension(:,:),intent(in)::S
  integer,intent(in)::i,j
  real(kind = 8)::left,right,top,bottom,spot
  spot = boundaries(s,i,j)
  top = boundaries(s,i,j + 1)
  bottom = boundaries(s,i, j - 1)
  right = boundaries(s,i + 1, j)
  left = boundaries(s,i - 1, j)

  delt = -0.5 * spot * ( top + bottom + left + right)

end function delt




real(kind = 8) function ToToNRG(s)
  implicit none
  real(kind = 8),intent(in),dimension(:,:)::s
  integer:: i,j
  real(kind = 8):: tot
  tot = 0.0
  do i = 1, size(s,2)
    do j = 1,size(s,2)
        tot = tot  + delt(s,i,j)
    end do
  end do
  ToToNRG = tot
end function ToToNRG

real(kind = 8) function magnetizm(s)
  implicit none
  real(kind = 8), dimension(:,:), intent(in)::s
  real(kind = 8) :: m
  integer:: i,j

  m = 0.0
  do i = 1, size(s,2)
    do j = 1, size(s,2)
      m = m + s(i,j)
    end do
  end do
  magnetizm = m/real(size(s,2)**2,kind = 8)
end function magnetizm


end program ising
