program fftwindow
use nr
use fftstuff
implicit none
real(kind = 8):: pi
integer::i,niquist
complex(kind =8), dimension(1024):: periodic,nonper,hanper,hannon,one,blacula,bnon
real(kind=8):: delta

niquist = 1024
delta =20.d0*3.141592d0/real(niquist,kind=8)
pi = 3.1415962
do i = 1,niquist
   periodic(i) = sin(2.d0 * pi * 10 * real(i,kind = 8)/real(niquist,kind = 8))
end do

open(12, file = "windowshit.txt",action='write')
open(13,file = 'blackmanstuff.txt',action = 'write')
one = 1.d0
nonper = 0
nonper(250:size(periodic)) = periodic(250:size(periodic))
do i = 1, niquist !index 0
   write(12,*)1.0*i/niquist,real(periodic(i))
end do
write(12,*)
write(12,*)
do i = 1, niquist !index 1
   write(12,*)1.0*i/niquist, real(nonper(i))
end do
write(12,*)
write(12,*)
call hanningwindow(periodic,hanper)
call hanningwindow(nonper,hannon)
call blackmann(nonper,blacula)
call four1(periodic,1)
call four1(nonper,1)
call fftflip(periodic)
call fftflip(nonper)

do i = 1, niquist !index 2
   write(12,*) - ts(1.d0,niquist,i),20.d0*log10(niquist**(-0.5)*abs(periodic(i))) 
end do
write(12,*)
write(12,*)
do i = 1, niquist !index 3
   write(12,*)ts(1.d0,niquist,i), 20.d0*log10(niquist**(-0.5)*abs(nonper(i)))
end do
write(12,*)
write(12,*)
do i = 1, niquist !index 4
   write(12,*)1.0*i/niquist,real(hanper(i))
end do
write(12,*)
write(12,*)
do i = 1, niquist !index 5
   write(12,*)1.0*i/niquist, real(hannon(i))
end do
write(12,*)
write(12,*)

do i = 1, niquist !index 6
   write(12,*)1.0*i/niquist, real(blacula(i))
end do
write(12,*)
write(12,*)

call four1(hanper,1)
call four1(hannon,1)
call four1(blacula,1)
call fftflip(blacula)
call fftflip(hanper)
call fftflip(hannon)
do i = 1, niquist !index 7
   write(12,*)  ts(1.d0,niquist,i),20.d0*log10(niquist**(-0.5)*abs(hanper(i)))
end do
write(12,*)
write(12,*)
do i = 1, niquist !index 8
   write(12,*)ts(1.d0,niquist,i), 20.d0 * log10(niquist**(-0.5)*abs(hannon(i)))
end do
write(12,*)
write(12,*)
do i = 1, niquist !index 9
   write(12,*)ts(1.d0,niquist,i), 20*log10(niquist**(-1.0)*abs(blacula(i)))
end do


close(12)
call blackmann(one,blacula)
do i = 1, niquist
   write(13,*) real(i),abs(blacula(i))
end do
write(13,*)
write(13,*)
call four1(blacula,1)
call fftflip(blacula)
do i = 1, niquist
   write(13,*) i, 10.d0 * log10(niquist**(1.0) * abs(blacula(i)))
end do
write(13,*)
write(13,*)
do i = 1, niquist
   write(13,*) -0.5 - 1.0 / niquist+ 1.0*i/niquist, 10.d0 * log10(niquist**(1.0) * abs(blacula(i)))
end do
write(13,*)
write(13,*)

   
close(13)
contains 

real(kind=8) function ts(dt,N,i)
  implicit none
  integer,intent(in)::N,i
  real(kind = 8)::dt

   ts= -(real(N,kind=8) + 2.d0*dt)/(2.d0) + real(i,kind = 8)/(dt) 
end function


end program fftwindow
