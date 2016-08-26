program filter
	!use fftstuff
	!use nr
	implicit none
	complex(kind = 4),dimension(3,4)::X
	complex(kind =4),dimension(4,3)::xt
	integer::i

	X(1,:) =(/1,2,3,9/)
	x(2,:) = (/8,5,1,2/)
	x(3,:) = (/9,8,7,2/)
	!x(4,:) = (/0,0,0,0/)
	
	call dft(x(:,1),1) !columns to start this bitch off
	call dft(x(:,2),1)
	call dft(x(:,3),1)
	call dft(x(:,4),1)

	do i = 1,3 !Print his bitch ass off. 
		print*, x(i,:)
	end do
	print*,
	print*,
	xt = transpose(x) !transpose that mother fucker
	!x(:,4) = 0
	call dft(xt(:,1),1) !Then do the columns again.
	call dft(xt(:,2),1)
	call dft(xt(:,3),1)
	!call dft(x(:,4),1)
	x = transpose(xt) !Then transpose that mother fucker AGAIN!
	

	do i = 1,3 !Print his bitch ass off. 
		print*, x(i,:)
	end do
	print*,
	print*,

contains
	!takes and transforms a rank 1 array vec and does the slow as balls DFT. Note it's intent in out, so your vector
	!Will be transformed. direction is either +1 for a forward transform, or -1 for an inverse transform. 
	!Everyt
	subroutine dft(vec,direction)
		implicit none
		complex(kind = 4),dimension(:),intent(inout):: vec
		integer,intent(in)::direction
		integer::i,j
		complex(kind = 4),dimension(size(vec),size(vec))::DFTM
		real::pi, n,norm
		n = 1.0 * size(vec)
		if (direction .ne. 1 .and. direction .ne. -1) then
			print*, "You dumb bastard, direction MUST be either +1 or -1. Integer. "
			print*,
			return
		else if(direction .eq. 1) then
			norm = 1.0
		else if (direction .eq. -1) then
			norm = 1.0/n
		end if 

		pi = 3.1415962
		
		do i =1,size(vec)
			do j = 1, size(vec)

				DFTM(i,j) = norm * exp(- real(direction)*2.0d0 * pi * (0,1) * (i-1)*(j-1)/n) 
				!This builds the DFT matrix. Note that (0,1) = sqrt(-1), norm is the normalization for 
				!inverse FFT, and direction is if you're doing an inverse or regular transform.
			end do
		end do
		
		vec = matmul(DFTM,vec)
		!print*, vec ,"IM VEC"
	end subroutine dft


	
end program filter