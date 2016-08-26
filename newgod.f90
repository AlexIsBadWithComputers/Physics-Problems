module ngf
	implicit none


contains
	!computes Godunov flux. This is a rewritten version of the codes available online at 
	!http://www.cfdbooks.com/cfdcodes/oned_euler_fluxes_v5.f90
	!Those suckers were used as guidance, but didn't work in their "as is formulation" so 
	!This is the Alex rewrite which made them work with gfortran.
	!also double precision. I also fixed some style things as well as removed some pointless
	!sections that didn't do what they needed to based on a direct copy paste. 
 
 function God(uL,uR)
		implicit none
		real(kind = 8)::God(3)
		real(kind = 8),intent(in):: uL(3), uR(3)
		real(kind = 8)::Flux(3)

		!constants for use in calculation
		real(kind = 8)::gamma,rhoL,rhoR,vL,vR,pL,pR !primitive variables
		real(kind =8)::al, aR !sound speed
		real(kind = 8)::Fp(3),Fm(3) !flux in both directions

		real(kind = 8):: gam,gam2,tol,pm1,mL,mR,vm,p(2),rm(2),r(2),rmI,amL,amR
		real(kind =8)::SmL,smR,Um2,Um3,um(3),pm2
		integer:: k,kmax

		!Constants
		gamma = 1.4
		tol = 10D-5
		gam2 = 0.5 * (gamma - 1.)/gamma

		!other variable declarations
		!stuff on the left
		rhoL=uL(1)
		vL = uL(2)/uL(1)
		pL = (gamma - 1.) * (uL(3) - 0.5 * rhoL * vL**2.)
		aL = sqrt(gamma * pR/rhoR)


		!stuff on right

		rhoR = uR(1)
		vR = uR(2)/uR(1)
		pR = (gamma - 1.) * (uR(3) - 0.5 * rhoR * vR**2)
		aR = sqrt(gamma * pL/rhoL)

		gam = (gamma + 1.)/(gamma - 1.)

		r(1) = rhoL
		r(2) = rhoR
		P(1) = pL
		P(2) = pR

		!Break into cases
		if (vL/aL .ge. 1.d0) then
			God = pf(uL)
			print*,"easy1"
			return
		else if (vR/aR .le. -1.d0) then
			God = pf(uR)
			print*, "easy"
			return 

		else 
			!Well, if the above didn't work (it probably didn't), we have to do it the hard way.

			pm1 =  (  (0.5*(vL-vR)*(gamma-1.)+aL+aR)/( &
           			  aL*pL**((1.-gamma)/gamma*0.5)   &
        	   + aR*pR**((1.-gamma)/gamma*0.5) )  )**(2.*gamma/(gamma-1.))

			!Now, solve that numerically 

			k = 0
			kmax = 100000

			do
				mL = massflux(rhoL,aL,pL,pm1)
				mR = massflux(rhoR,aR,pR,pm1)

				pm2 = ( mL * pR + mR * pL - mL * mR * (vR - vL))/(mL + mR)

				k = k + 1

				if( abs( pm2 - pm1) .lt. tol) exit
				if ( k .gt. kmax) then
					print*, "You broke me. You bastard., line 85 of God.f90"
					stop 
				end if 

				pm1 = pm2
			end do

			mL = massflux(rhoL,aL,pL,pm2)
			mR = massflux(rhoR,aR,pR,pm2)
			vm = (mL*vL+mR*vR - (pR-pL))/(mL+mR)

			!density in betweenies

			do k = 1,2
				if (pm2/p(k) .ge.1.d0) then
					rm(k) = r(k) *( 1.d0 + gam * pm2/p(k))/(gam + pm2/p(k))
				else
					rm(k) = r(k) * (pm2/p(k)) **( 1./gamma)
				end if  
			end do 

			!Where do the waves make contact?

			if(vm .ge. 0.d0) then
				rmI = rm(1)
			else
				rmI = rm(2)
			end if  

			!Find wave speeds
			amL = sqrt(gamma * pm2/rm(1))
			amR = sqrt(gamma * pm2/rm(2))
			smL = vm - amL
			smR = vm + amR

			!sonic case

			if(SmL .le. 0.d0 .and. smR .ge. 0.d0) then
				um2 = rmi * vm
				um3 = Pm2/(gamma - 1.) + 0.5 * rmI * vm**2
			else if( SmL .gt. 0.d0 .and. vR + aR .gt. 0.d0) then
				call sonic(rmI,Um2,Um3,vR,aR,PR,vm,amR,vR+aR,SmR)
  			 endif

  			 !compute flux

  			 um(1) = rmI
  			 um(2) = um2
  			 um(3) = um3
  			 God = pf(um)

		end if

	end function

	
	
!From CFD is fun (same source as above.)
!again this is the alex version with a few style changes and "whoopsie" fixes
!as compared to the online version.

function VanLeer(uL,uR)
 implicit none
 real(kind = 8) :: uL(3), uR(3) !  Input (conservative variables rho*[1, v, E])
 real(kind = 8) :: VanLeer(3)   ! Output (numerical flux across L and R states)
!Local constants
 real(kind = 8) :: gamma                        ! Ratio of specific heat.
!Local variables
 real(kind = 8) :: rhoL, rhoR, vL, vR, pL, pR   ! Primitive variables.
 real(kind = 8) :: aL, aR, ML, MR               ! Speeds of sound and Mach numbers.
 real(kind = 8) :: Fp(3), Fm(3)                 ! F+ and F-
!Constants.
     gamma = 1.4
     

!Primitive and other variables.
!  Left state
    rhoL = uL(1)
      vL = uL(2)/uL(1)
      pL = (gamma-1.0d0)*( uL(3) - 0.5d0*rhoL*vL*vL )
      aL = sqrt(gamma*pL/rhoL)
      ML = vL/aL !left mach number

!  Right state
    rhoR = uR(1)
      vR = uR(2)/uR(1)
      pR = (gamma-1.0d0)*( uR(3) - 0.5d0*rhoR*vR*vR )
      aR = sqrt(gamma*pR/rhoR)
      MR = vR/aR   !right mach numer

!Positive Part of Flux evaluated in the left cell.
 Fp(1) =   0.25d0*rhoL*aL*(ML+1.0d0)*(ML+1.0d0)
 Fp(2) =   Fp(1)*2.0d0*aL*(1.0d0+0.5d0*(gamma-1.0d0)*ML)/gamma
 Fp(3) =   Fp(1)*2.0d0*aL*aL*(1.0d0+0.5d0*(gamma-1.0d0)*ML)**2/(gamma*gamma-1.0d0)

!Negative Part of Flux evaluated in the right cell.
 Fm(1) = - 0.25d0*rhoR*aR*(MR-1.0d0)*(MR-1.0d0)
 Fm(2) =   Fm(1)*2.0d0*aR*(-1.0d0+0.5d0*(gamma-1.0d0)*MR)/gamma
 Fm(3) =   Fm(1)*2.0d0*aR*aR*(1.0d0-0.5d0*(gamma-1.0d0)*MR)**2/(gamma*gamma-1.0d0)

!Compute the flux: Fp(uL)+Fm(uR).
   VanLeer = Fp + Fm

 end function VanLeer





	 function pf(u)
 		real(kind = 8) :: u(3)             !  Input (conservative variables [rho, rho*v, rho*E])
 		real(kind = 8) :: pf(3) ! Output (physical flux of the Euler equations)
		!Local variables
		real(kind = 8) :: density, velocity, pressure, enthalpy, gamma

		!Define and compute some quantities.
 	    gamma = 1.4
	    density = u(1)
  		velocity = u(2)/u(1)
		pressure = (gamma-1.0)*( u(3) - 0.5*density*velocity*velocity )
		enthalpy = u(3) + pressure

		!Evaluate the physical flux (mass, momentum, and energy fluxes).
		pf(1) =           density * velocity
 		pf(2) = (density*velocity)* velocity + pressure
 		pf(3) =          enthalpy * velocity

 end function pf


function massflux(r,c,pQ,pm)
 	 real(kind = 8) :: r,c,pQ,pm !  Input
 	 real(kind = 8) :: massflux  ! Output
	 !Local variables
  	 real(kind = 8) :: gamma,gam1,gam2,eps
   
     gamma = 1.4
     eps = 1.0e-10

     gam1=0.5d0*(gamma+1.0d0)/gamma; gam2=0.5d0*(gamma-1.0d0)/gamma;
     if (pm/pQ >= 1.0d0-eps) then ! eps to avoid zero-division
     	 massflux=(r*c)*sqrt( 1.0d0+gam1*(pm/pQ-1.0d0) )
     else
     	 massflux=(r*c)*gam2*(1.0d0-pm/pQ)/( 1.0d0-(pm/pQ)**(gam2) )
     endif
  end function massflux

	subroutine sonic(US1,US2,US3, u1,c1,P1,u2,c2,a1,a2)
  		real(kind = 8) :: u1,c1,P1,u2,c2,a1,a2 ! Input
  		real(kind = 8) :: US1,US2,US3          ! Output
		!Local variables
  		real(kind = 8)              :: us,cs,Ps,rs,R1,R2    
  		real(kind = 8) :: gamma
  		gamma = 1.4
  		

   		R1 =  a2/(a2-a1)
    	R2 = -a1/(a2-a1)
   		us = R1*u1+R2*u2
   		cs = R1*c1+R2*c2
   		Ps = (cs/c1)**(2.0d0*gamma/(gamma-1))*P1
 	    rs = gamma*Ps/(cs*cs)

  		US1 = rs
 		US2 = rs*us
  	    US3 = Ps/(gamma-1.0d0)+0.5d0*rs*us*us

  end subroutine sonic


!This is my accidental lax-wendrov. I wrote it by accident and I figured why not use it.


subroutine LW(gp,uf,u,c,p,totaltime)
		implicit none
		integer::i,j,n,k
		real(kind=8),intent(in)::totaltime
		integer,intent(in)::gp! number of grid points
		real(kind = 8):: t,gamma,tmax,length,r,pr,e,press,dt,dx
		real(kind = 8),dimension(gp + 2),intent(out):: p,c,u !u is speed, c is sound speed
		real(kind = 8),dimension(gp+ 2,3), intent(out)::uf
		real(kind =8),dimension(gp + 2,3):: newuf
		real(kind = 8),dimension(gp + 2,3)::F   !fluxes       
		real(kind = 8), parameter :: eps=1E-15, CFL = 0.5



		!initial conditions/known things
		
		gamma = 1.4d0
		!fill from the left
		do i = 1,nint(1.0*gp/2) !boundary condition are farther cause it big
			uf(i,1) = 10.**5
			p(i) = 1.0
			uf(i,2) = uf(i,1) * u(i)
			u(i) = 0.0
			c(i) = sqrt(gamma * p(i)/uf(i,1))
			uf(i,3) = p(i)/(gamma-1.)
		end do
			
		do i = nint(1.0*gp/2)+1,gp +  2
			uf(i,1) = 1.25 * 10 ** 4
			p(i) = 0.1
			u(i) = 0.0
			uf(i,2) = uf(i,1) * u(i)
			uf(i,3) =p(i)/(gamma-1.) 
			c(i) = sqrt(gamma * p(i)/uf(i,1))

		end do 
		length = 1.0
		dx =.01! length/(gp -1.)
		t = 0. 
		F = 0.
		newuf = 0.

		
		do while(t .lt. totaltime) 
			dt = CFL * dx/maxval(c(:))

			t = t + dt
					
			do j = 1,gp+1
				r = uf(j,1)
				pr = uf(j,2) !pr is rho r, the second term in the equations
				e = uf(j,3)
				press = (gamma - 1.) * (e - 0.5 * pr * pr / r  )
				F(j,1) = pr 
				F(j,2) = pr * pr / r + press
				F(j,3) = pr / r * (e + press) 
			
			end do

			!Now we need to find the half steps
			do j = 1,gp+1
				do i = 1,3
					newuf(j,i) =  0.5 * (uf(j + 1, i) + uf(j, i)) - 0.5 * dt/dx * (F(j + 1,i) - F(j,i))
				end do
			end do
			call reflectivebc(newuf,gp+2)
			

			!Now, we need half step of flux for the next part
			do j = 1,gp+1
				r = newuf(j,1)
				pr = newuf(j,2) !pr is rho r, the second term in the equations
				e = newuf(j,3)
				press = (gamma - 1.) * (e - 0.5 * pr * pr / r  )

				F(j,1) = pr 
				F(j,2) = pr * pr / r + press
				F(j,3) = pr / r * (e + press) 
				
			end do
				
			!Now we need to update using those half steps as well
			do j = 2,gp
				do i = 1,3
					newuf(j,i) =uf(j,i) - 0.5 * dt/dx * (F(j,i) - F(j-1,i))
				end do
			end do
			!Now we need to update my man uf

			do j =2,gp
				do i = 1,3
					uf(j,i) = newuf(j,i)
				end do
		!	print*, uf(j,1), uf(j,2), uf(j,3),k,j
			end do

			do j = 2,gp
				u(j) = uf(j,2)/uf(j,1)    !calculate your other shit.
				p(j) = press
				c(j) = sqrt(gamma * abs(p(j))/uf(j,1))
			end do	
			print*, maxval(c(:)),dt,"DFD"
		end do		
		
	end subroutine LW

	
	!This is for the Lax wenrov
	subroutine reflectivebc(U,n)
		implicit none
		integer,intent(in)::n
		real(kind = 8),dimension(n,3),intent(inout) :: U

		U(1,1) = U(2,1)
		U(1,2) = - U(2,2)
		U(1,3) = U(2,3)
		U(n-1,1) = U(n-2,1) 
		U(n-1,2) = -U(n-2,2) 
		U(n-1,3) = U(n-2,3) 
	end subroutine




end module