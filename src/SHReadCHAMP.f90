subroutine SHReadCHAMP(filename, cilm, lmax, gm, r0_pot, error)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This program will read the file format of spherical harmonic
!	coefficients used by the CHAMP/GRACE group into standard arrays for
!	of Cilm and the error.
!
!	Calling Parameters
!		IN
!			filename: 	The name of the file.
!		OUT
!			lmax:		Maximum spherical harmonic degree.
!			cilm:		An array of the spherical harmonic coefficients.
!			gm:		GM
!			R0_pot		Reference radius of gravity field
!		OPTIONAL
!			error:	An array containing the error coefficients.
!
!	Written by Mark Wieczorek (September 2005)
!
!	Copyright (c) 2005, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none
	
	character(*), intent(in) ::	filename
	integer, intent(out) ::	lmax
	real*8, intent(out) ::	cilm(:,:,:), gm, r0_pot
	real*8, intent(out), optional :: error(:,:,:)
	integer ::	temp_max
		
	integer:: l, m, stat
	character :: c*6
	real*8 ::	clm, slm, clmsig, slmsig
	
	cilm=0.0d0
	
	if (size(cilm(:,1,1)) < 2) then
		print*, "Error --- SHReadCHAMP"
		print*, "CILM must be dimensioned (2, *, *)."
		print*, "Input array is dimensioned ", size(cilm(:,1,1)), size(cilm(1,:,1)), size(cilm(1,1,:))
		stop
	endif
	
	if (present(error)) then
		if (size(error(:,1,1)) < 2) then
			print*, "Error --- SHReadCHAMP"
			print*, "ERROR must be dimensioned (2, *, *)."
			print*, "Input array is dimensioned ", size(error(:,1,1)), size(error(1,:,1)), size(error(1,1,:))
			stop
		endif
		
		error=0.0d0
	endif
	
	
	open(13,file=filename)
	
	lmax = 0
	temp_max = 0
	gm = 0.0d0
	r0_pot = 0.0d0
	
	do 
		read(13,"(a6)", advance="no", iostat=stat) c
		if(stat /=0) exit
		
		if (c=="EARTH " .or. c(1:3) =="GGM") then
			read(13,*) gm, r0_pot
		elseif(c=="SHM   ") then 
			read(13,*) lmax
			if (size(cilm(1,:,1)) < lmax+1 .or. size(cilm(1,1,:)) < lmax+1) then
				print*, "Error ---SHReadCHAMP"
				print*, "CILM must be dimensioned (2, LMAX+1, LMAX+1) where LMAX is ", lmax
				print*, "Input array is dimensioned ", size(cilm(1,:,1)),  size(cilm(1,1,:)) 
				stop
			elseif (present(error)) then
				if (size(error(1,:,1)) < lmax+1 .or. size(error(1,1,:)) < lmax+1) then
					print*, "Error ---SHReadCHAMP"
					print*, "ERROR must be dimensioned (2, LMAX+1, LMAX+1) where LMAX is ", lmax
					print*, "Input array is dimensioned ", size(error(1,:,1)),  size(error(1,1,:)) 
					stop
				endif
			endif
		elseif(c=="GRCOF2" .or. c=="CALSDV" .or. c=="gfc") then
			if (present(error)) then
				read(13,*) l, m, clm, slm, clmsig, slmsig
				cilm(1,l+1,m+1) = clm
				cilm(2,l+1,m+1) = slm
				error(1,l+1,m+1) = clmsig
				error(2,l+1,m+1) = slmsig
			else
				read(13,*) l, m, clm, slm
				cilm(1,l+1,m+1) = clm
				cilm(2,l+1,m+1) = slm
			endif
			
			if (l > temp_max) temp_max = l
		else
			read(13,*)
			print*, "SHReadCHAMP --- Ignoring line ", c
		endif
	enddo
	
	if (lmax == 0) lmax = temp_max

	close(13)
	
end subroutine SHReadCHAMP
