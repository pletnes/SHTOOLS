subroutine Hilm(cilm, ba, gridglq, lmax, nmax, mass, r0, rho, w, plx, zero, filter_type, filter_deg)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This routine will compute the next estimate of Moho coefficients 
!	given an initial estimate of the moho in a gridded data file. This
!	is simply Equation 18 in Wieczorek and Phillips (1998). Note that the
!	degree-0 topography must be included in the gridded relief. Note that 
!	the array plx is optional, and should not be precomputed when memory 
!	is an issue (i.e., lmax>360).
!
!	Calling Parameters:
!		IN
!			ba:		Bouguer Anomaly spherical harmonic coefficients.
!			gridglq:	Initial estimate of Moho refief, gridded according
!					to a call to MakeGridGLQ.
!			lmax:		Maxmimum spherical harmonic degree to compute.
!			nmax:		Order of potential coefficient expansion.
!			mass:		Mass of planet.	
!			R0:		The radius that the coefficients are referenced
!					to. 		
!			rho:		Density contrast between the mantle and crust.
!			w:		Gauss-Legendre points used in integrations
!					(determined from a call to PreCompute).
!		OUT
!			cilm:		Estimate of Moho relief spherical harmonic coefficients, 
!					with dimensions (2, lmax+1, lmax+1).
!			
!		OPTIONAL		
!			filter_deg:	If specified, each interation will be filtered according to
!					equation 18 in Wieczorek and Phillips (1998), where the value
!					of filter_deg corresponds to the spherical harmonic degree where 
!					the filter is 1/2.
!			filter_type:	If filter_deg is specified, this must be as well. 
!					A value of (1) corresponds to the minimum amplitude filter in 
!					Wieczorek and Phillips (1998), whereas as value of (2) corresponds
!					to a minimum curvature filter.
!			plx:		Input array of Associated Legendre Polnomials computed
!					at the Gauss points (determined from a call to
!					PreCompute). If this is not included, then the optional
!					array zero MUST be inlcuded.
!			zero		Array of dimension lmax+1 that contains the latitudinal
!					gridpoints used in the Gauss-Legendre quadrature integration
!					scheme. Only needed if plx is not included.
!
!	All units assumed to be SI.
!
!	Dependencies:		NGLQSH, SHExpandGLQ, wl, wl_curv
!
!	Written by Mark Wieczorek 2003
!	September 3, 2005. Modifed so that the array plx is now optional.
!
!	Copyright (c) 2005, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	use SHTOOLS, only: NGLQSH, SHExpandGLQ, Wl, WlCurv

	implicit none
		
	real*8, intent(out) :: 			cilm(:,:,:)
	real*8, intent(in) ::			ba(:,:,:), gridglq(:,:), mass, r0, rho, w(:)
	real*8, intent(in), optional ::		plx(:,:), zero(:)
	integer, intent(in) ::			lmax, nmax
	integer, intent(in), optional :: 	filter_type, filter_deg
	real*8 ::				prod, pi, d, depth, filter(lmax+1)
	real*8, allocatable ::			cilmn(:, :, :), grid2(:,:)
	integer  ::				j, l, n, nlong, nlat, astat(2)


	if (size(cilm(:,1,1)) < 2 .or. size(cilm(1,:,1)) < lmax+1 .or. size(cilm(1,1,:)) < lmax+1) then
		print*, "Error --- HilmGLQ"
		print*, "CILM must be dimensioned as (2, LMAX+1, LMAX+1) where LMAX is ", lmax
		print*, "Input dimension is ", size(cilm(:,1,1)), size(cilm(1,:,1)), size(cilm(1,1,:))
		stop
	elseif (size(ba(:,1,1)) < 2 .or. size(ba(1,:,1)) < lmax+1 .or. size(ba(1,1,:)) < lmax+1) then
		print*, "Error --- HilmGLQ"
		print*, "BA must be dimensioned as (2, LMAX+1, LMAX+1) where LMAX is ", lmax
		print*, "Input dimension is ", size(ba(:,1,1)), size(ba(1,:,1)), size(ba(1,1,:))
		stop
	elseif(size(gridglq(1,:)) < 2*lmax+1 .or. size(gridglq(:,1)) <  lmax +1) then
		print*, "Error --- HilmGLQ"
		print*, "GRIDGLQ must be dimensioned as (LMAX+1, 2*LMAX+1) where LMAX is ", lmax
		print*, "Input dimension is ", size(gridglq(1,:)), size(gridglq(:,1))
		stop
	elseif(size(w) < lmax+1) then
		print*, "Error --- HilmGLQ"
		print*, "W must be dimensioned as (LMAX+1) where LMAX is ", lmax
		print*, "Input dimension is ", size(w)
		stop
	endif
	
	if (present(zero)) then
		if (size(zero) < lmax + 1) then
			print*, "Error --- HilmGLQ"
			print*, "ZERO must be dimensioned as (LMAX+1) where LMAX is ", lmax
			print*, "Input dimension is ", size(zero)
			stop
		endif
	endif
	
	if (present(plx)) then
		if (size(plx(:,1)) < lmax+1 .or. size(plx(1,:)) < (lmax+1)*(lmax+2)/2) then
			print*, "Error --- HilmGLQ"
			print*, "PLX must be dimensioned as (LMAX+1, (LMAX+1)*(LMAX+2)/2) where LMAX is ", lmax
			print*, "Input dimension is ", size(plx(:,1)), size(plx(1,:))
			stop
		endif
	endif
	
	if(present(filter_type)) then
		if (filter_type /=0 .and. filter_type /=1 .and. filter_type /=2) then
			print*, "Error --- HilmGLQ"
			print*, "FILTER_TYPE must be either 0 (none), 1 (minimum amplitude), or 2 (minimum curvature)"
			print*, "Input value is ", filter_type
			stop
		endif
	endif
	
	if (present(plx)) then
		continue
	else
		if (present(zero)) then
			continue
		else
			print*, "Error --- Hilm"
			print*, "If the optional array PLX is not included, it is necessary to include the optional array ZERO."
			stop
		endif
	endif
	
	allocate(cilmn(2, lmax+1, lmax+1), stat = astat(1))
	allocate(grid2(lmax+1,2*lmax+1), stat = astat(2))
	if (astat(1) /= 0 .or. astat(2) /= 0) then
		print*, "Error --- HilmGLQ"
		print*, "Problem allocating arrays CILMN and GRID2", astat(1), astat(2)
		stop
	endif

	cilm = 0.0d0
	cilmn = 0.0d0
	pi = acos(-1.0d0)
	nlat = NGLQSH(lmax)
	nlong = 2*lmax+1
	grid2(1:nlat,1:nlong) = gridglq(1:nlat,1:nlong)
	
	if (present(plx)) then
		call SHExpandGLQ(cilmn, lmax, grid2(1:nlat,1:nlong), w(1:lmax+1), plx = plx(1:lmax+1, 1:(lmax+1)*(lmax+2)/2), &
			norm = 1, csphase = 1)
	else 
		call SHExpandGLQ(cilmn, lmax, grid2(1:nlat,1:nlong), w(1:lmax+1), zero = zero(1:lmax+1), norm = 1, csphase = 1)
	endif

	d = cilmn(1,1,1)
	depth = r0-d
	print*, "Average depth of Moho (km) = ", depth/1.d3
	
	cilm(1,1,1) = d
	grid2(1:nlat,1:nlong) = grid2(1:nlat,1:nlong) - d
	
	filter = 1.0d0
	if (present(filter_type) .and. present(filter_deg)) then
		do l=1, lmax
			if (filter_type ==1) then
				filter(l+1) = Wl(l, filter_deg, r0, d)
			elseif(filter_type==2) then
				filter(l+1) = WlCurv(l, filter_deg, r0, d)
			endif
		enddo
	endif
	
	do l=1, lmax
		cilm(:,l+1,1:l+1) = filter(l+1)*ba(:,l+1,1:l+1) * mass * dble(2*l+1) * ( (r0/d)**l) / &
			(4.0d0*pi*rho*d**2)
	enddo
		
	do n=2, nmax
	
		if (present(plx)) then
			call SHExpandGLQ(cilmn, lmax, (grid2(1:nlat,1:nlong)/d)**n, w(1:lmax+1), &
				plx = plx(1:lmax+1, 1:(lmax+1)*(lmax+2)/2), norm = 1, csphase = 1)
		else
			call SHExpandGLQ(cilmn, lmax, (grid2(1:nlat,1:nlong)/d)**n, w(1:lmax+1), zero = zero(1:lmax+1), &
				norm = 1, csphase = 1)
		endif
		
		do l = 1, lmax
			prod = 1.0d0
			do j=1, n
				prod = prod * dble(l+4-j)
			enddo
			prod = d * prod/( dble(l+3) * dble(fact(n)) )
			
			cilm(:,l+1,1:l+1) = cilm(:,l+1,1:l+1) - filter(l+1)*cilmn(:,l+1,1:l+1)*prod
		enddo
	enddo
	
	deallocate(cilmn)
	deallocate(grid2)	
	
	contains
	
		function fact(i)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!
		!	This function computes the factorial of an integer
		!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			implicit none
			integer ::	i, j, fact
	
			if (i == 0) then
				fact = 1	
			elseif (i .lt. 0) then
				print*, "Argument to FACT must be positive"
				stop
			else
				fact = 1
				do j = 1, i
					fact = fact * j
				enddo
			endif
			
		end function fact
	
end subroutine Hilm

