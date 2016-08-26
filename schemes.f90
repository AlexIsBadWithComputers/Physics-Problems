program schemes
!	use other
	use godunov_flux
	implicit none
	real(kind = 8),dimension(1002) :: p, u,c
	real(kind = 8), dimension(1002,3) :: uf
	real(kind = 8):: u2(1002,3)
	integer :: n , i 
	n = 1000
		open(12,file = "output.txt",action = "write")
	call LW(n,uf,u,c,p,10000.0d0)
	print*,"done"
	
	
	call useRiemann(n,U2)
		do i = 3,size(U2(:,1))
			!print*,i  , u2(i,1), (1.4 - 1.) * (u2(i,3) - 0.5 * u2(i,2) * u2(i,2) / u2(i,1)), u2(i,3)/u2(i,1)

		write(12,*) 0.1*i   , u2(i,1), (1.4 - 1.) * (u2(i,3) - 0.5 * u2(i,2) * u2(i,2) / u2(i,1)),&
		 u2(i,3)/(u2(i,1)) - u2(i,2) * u2(i,2) / u2(i,1) / u2(i,1)/2.,u2(i,2)/u2(i,1)
	end do

	write(12,*) 
	write(12,*)
	print*, "van start"
	do i = 2,size(Uf(:,1))
			
		write(12,*) 1.0*i/10.   , uf(i,1), (1.4 - 1.) * (uf(i,3) - 0.5 * uf(i,2) * uf(i,2) / uf(i,1)),& 
		uf(i,3)/(uf(i,1)) - uf(i,2) * uf(i,2) / uf(i,1) / uf(i,1)/2. ,uf(i,2)/uf(i,1)
	end do

	close(12)
contains
	


	subroutine useRiemann(gp,U)
		implicit none
		integer,intent(in)::gp
		real(kind = 8), intent(out),dimension(gp + 2,3) :: U
		
		real(kind = 8),dimension(gp + 2,3)::F,un
		integer :: j,i,n,m
		real(kind = 8)::t,dt,dx,CFL,p,gamma,dumF(3),c(gp),press,v(gp)
		CFL = 0.5
		dx = .12
		t = 0.0
		gamma = 1.4
		n = 0
		m = 0
		Un = 0
		do i = 1,nint(1.0*gp/2) !boundary condition are farther cause it big
			u(i,1) = 100000.
			p = 1.0
			u(i,2) = 0.0!uf(i,1) * v
			u(i,3) = p/(gamma-1.)
		end do
			
		do i = nint(1.0*gp/2)+1,gp +  2
			u(i,1) = 1.25E4
			p = 0.1
			u(i,2) =0.! uf(i,1) * v
			u(i,3) =p/(gamma-1.) 
		end do

		do j = 1,gp-2
				press = (gamma - 1.) *( u(j,3) - 0.5 * u(j,2) * u(j,2) / u(j,1))
				c(j) = sqrt(gamma *press/u(j,1))
		end do
	
		do while (t .lt. 5000.)
			dt = cfl * dx / maxval(c)

			t = t + dt
			
			
			do j = 2, gp
				do i = 1,3
					U(j,i) = 0.5*(U(j+1,i) + U(j-1,i))
				end do 
			end do 
			
			do j = 2,gp
				dumF =  vanleer(U(j-1,:),U(j,:))
				F(j,1) = dumF(1)
				F(j,2) = dumF(2)
				F(j,3) = dumF(3)
			end do
			
			do j = 2, gp
				do i = 1,3
					!F(j,i) = 0.5*(F(j+1,i) + F(j,i))
				end do 
			end do 

			do j = 2,gp
				do i = 1,3
					U(j,i) = U(j,i) - dt/dx * (F(j,i) - F(j-1,i))
				end do	
				!print*, u(j,1)
			end do
			
			do j = 2,gp
				press = (gamma - 1.) *( u(j,3) - 0.5 * u(j,2) * u(j,2) / u(j,1))
				v(i) = u(i,2)/u(i,1)
				c(j) = sqrt(abs(gamma *press/u(j,1)))
				!print*, c(j)
			end do
		
		end do
		print*, "Done motherfucker."
			
	end subroutine





			
	
end program