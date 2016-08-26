program binary
	use nr
	use fftstuff
	implicit none
	real(kind = 4), allocatable,dimension(:):: x,y
	real(kind = 4),dimension(:),pointer :: px,py
	real(kind = 4), dimension(:), allocatable::px2,pxy
	real:: T,prb,epy,m
	integer:: jmax,lines,i
	real(kind = 4):: prob,a,b
	!allocate(px(10),py(10)
	call filelength("binary.txt",lines)	
	!print*, lines, "Im here"
	open(12,file = 'binary.txt', action = 'read')
	allocate(x(lines),y(lines))
	do i= 1, lines
		read(12,*) A,B
		x(i) = a
		y(i) = b
	end do
	
	T = maxval(x)-minval(x)
	print*,T, minval(x), maxval(x), 1.0 * lines/2./T

	
	call fasper(x,y,.5,4.,px,py,jmax,prob)
	!call fftflip2(px)
	!call fftflip2(py)
	
	
	print*, jmax, prob, maxval(px),lines - size(px)
	
	open(12,file = "binaryfasper.txt",action = 'write')
print*, px(jmax),jmax ,"yo yo yo"
	!x = x(size(x):1:-1)
	
	do i = 1,size(py)
		write(12,*) px(i), 10*log10(py(i))

	end do
	write(12,*) 
	write(12,*)
	do i = 1, size(py)
		write(12,*)px(i), (py(i))
	end do
	write(12,*) 
	write(12,*)
	
	m = 0.5 * 2.0 * size(py)
	do i = 1, size(py)
		epy = exp(-py(i))
		prb = epy*m
		if (prb > 0.00001) then
			 prb = 1.0 - (1. - epy)**m
		end if
		if (prb < 0.01) then
			print*, prb, px(i),i,py(i),log(py(i))
		end if 
		if (prb .gt. 0.30 .and. prb .lt. 0.8) then 
			print*,
			print*, "Bottom Acceptance Rate"
			print*, prb, px(i),i, py(i), log(py(i))
			print*,
		end if
	end do 










end program binary