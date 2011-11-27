! Copyright 2011, Paul Anton Letnes

! Functions and subroutines that do not need wrapping:
!   Random.f95: Replace by 'import random' or 'import numpy.random' in python
!   FFTW3.f95:  Seems to be 'internal use only', no need to expose it

module pyshtools
    implicit none
    complex(8), parameter :: imunit = (0.0_8, 1.0_8)
    private :: imunit

contains

! Various support functions
function PlmIndex(l, m)
    use shtools, only: PlmIndex_f => PlmIndex
    implicit none
    integer :: PlmIndex
    integer, intent(in) :: l, m

    ! TODO check fortran -> python index issue.
    ! Fortran: 1-based. Python: 0-based
    PlmIndex = PlmIndex_f(l, m) - 1
end function PlmIndex

function YilmIndex(i, l, m)
    use shtools, only: YilmIndex_f => YilmIndex
    implicit none
    integer :: YilmIndex
    integer, intent(in) :: i, l, m

    ! TODO check fortran -> python index issue.
    ! Fortran: 1-based. Python: 0-based
    YilmIndex = YilmIndex_f(i, l, m) - 1
end function YilmIndex

! Legendre functions
subroutine PLegendre(p, lmax, z)
    use shtools, only: PLegendre_f => PLegendre
    implicit none
    integer, intent(in) :: lmax
    real(8), intent(out) :: p(lmax + 1)
    real(8), intent(in) :: z
    !f2py depend(lmax) p

    call PLegendre_f(p, lmax, z)
end subroutine PLegendre


subroutine PlBar(p, lmax, z)
    use shtools, only: PlBar_f => PlBar
    implicit none
    integer, intent(in) :: lmax
    real(8), intent(out) :: p(lmax + 1)
    real(8), intent(in) :: z
    !f2py depend(lmax) p

    call PlBar_f(p, lmax, z)
end subroutine PlBar

subroutine PlSchmidt(p, lmax, z)
    use shtools, only: PlSchmidt_f => PlSchmidt
    implicit none
    integer, intent(in) :: lmax
    real(8), intent(out) :: p(lmax + 1)
    real(8), intent(in) :: z
    !f2py depend(lmax) p
    
    call PlSchmidt_f(p, lmax, z)
end subroutine PlSchmidt

subroutine PlON(p, lmax, z)
    use shtools, only: PlON_f => PlON
    implicit none
    integer, intent(in) :: lmax
    real(8), intent(out) :: p(lmax + 1)
    real(8), intent(in) :: z
    !f2py depend(lmax) p
    
    call PlON_f(p, lmax, z)
end subroutine PlON

! Associated Legendre functions
subroutine PlmBar(p, lmax, z, csphase, cnorm)
    use shtools, only: PlmBar_f => PlmBar
    implicit none
    integer, intent(in) :: csphase, cnorm
    integer, intent(in) :: lmax
    real(8), intent(out) :: p((lmax + 1)*(lmax + 2) / 2)
    real(8), intent(in) :: z
    !f2py depend(lmax) p
    !f2py integer, optional :: csphase=1
    !f2py integer, optional :: cnorm=0

    call PlmBar_f(p, lmax, z, csphase, cnorm)
end subroutine PlmBar

subroutine PlmSchmidt(p, lmax, z, csphase, cnorm)
    use shtools, only: PlmSchmidt_f => PlmSchmidt
    implicit none
    integer, intent(in) :: csphase, cnorm
    integer, intent(in) :: lmax
    real(8), intent(out) :: p((lmax + 1)*(lmax + 2) / 2)
    real(8), intent(in) :: z
    !f2py depend(lmax) p
    !f2py integer, optional :: csphase=1
    !f2py integer, optional :: cnorm=0

    call PlmSchmidt_f(p, lmax, z, csphase, cnorm)
end subroutine PlmSchmidt

subroutine PlmON(p, lmax, z, csphase, cnorm)
    use shtools, only: PlmON_f => PlmON
    implicit none
    integer, intent(in) :: csphase, cnorm
    integer, intent(in) :: lmax
    real(8), intent(out) :: p((lmax + 1)*(lmax + 2) / 2)
    real(8), intent(in) :: z
    !f2py depend(lmax) p
    !f2py integer, optional :: csphase=1
    !f2py integer, optional :: cnorm=0

    call PlmON_f(p, lmax, z, csphase, cnorm)
end subroutine PlmOn

! Subroutines written by me, for me)
!subroutine ylm(y, lmax, theta, phi)
    !! Calculate the spherical harmonics for all l <= lmax and for
    !! corresponding m values.
    !use shtools, only: PlmSchmidt, PlmIndex
    !implicit none
    !integer, intent(in) :: lmax
    !complex(8), dimension((lmax + 1) * (lmax + 1)), intent(out) :: y
    !real(8), intent(in) :: theta, phi
    !!f2py depend(lmax) y
    !! Local variables
    !real(8), dimension((lmax+1) * (lmax+2) / 2) :: p
    !real(8) :: invsqrt4pi
    !complex(8), dimension(lmax) :: mfacs
    !integer :: ints(lmax)
    !integer :: l, m, i
    !integer :: signfac = 1
    !invsqrt4pi = 1.0_8 / sqrt(4.0_8 * 3.1415926535897932385_8)

    !forall (i = 1 : size(ints))
        !ints(i) = i
    !end forall

    !! Fill the mfacs array with \exp(i m \phi) phase factors
    !mfacs = cos(ints(: lmax) * phi) + imunit * sin(ints(: lmax) * phi)

    !! Fill up the array p with the values of the assoc. Legendre func.
    !! cnorm=1: the complex normalization, see SHTOOLS documentation
    !call PlmSchmidt(p, lmax, cos(theta), cnorm=1)

    !do l = 0, lmax
        !! First, fill the m = 0 spot, to get a symmetric l-loop
        !y(YlmIndex(l, 0)) = invsqrt4pi * sqrt(2.0_8 * l + 1.0_8) &
            !* p(PlmIndex(l, 0))
        !! signfac stores the sign of the factor (-1)**m
        !signfac = 1
        !do m = 1, l
            !signfac = signfac * (-1)
            !i = YlmIndex(l, m) 
            !y(i) = signfac * invsqrt4pi &
                !* sqrt(2.0_8 * l + 1.0_8) * p(PlmIndex(l, m)) * mfacs(m)
            !y(YlmIndex(l, -m)) = signfac * conjg(y(i))
        !end do
    !end do
!end subroutine ylm

!pure function ylmindex(l, m)
    !! Calculate the index of a given Ylm multipole in the array.
    !! NB: Internal use only!
    !implicit none
    !integer :: ylmindex
    !integer, intent(in) :: l, m

    !ylmindex = l**2 + l + m + 1
!end function ylmindex

end module

