MODULE void_parameters
    implicit none

    integer, parameter :: nx = 12           !number of variables in the system of equations
    !SPACE & TIME DOMAIN RESOLUTION:
    integer, parameter :: nr = 200          !number of radial grid points
    integer, parameter :: nt = 100000          !number of time grid points, might make this variable later. 
    !Cosmological Model Parameters:
    double precision, parameter :: omega_matter = 0.3
    double precision, parameter :: omega_lambda = 1.0 - omega_matter
    double precision, parameter :: H0 = 67.810d0 !
    !
    !double precision :: age
    double precision, parameter :: pi = 4d0*datan(1.0d0)
    !non-dimensioanlisations/units
    double precision, parameter :: mu = 1.989d45            ! mass in 10^15 M_{\odot}
    double precision, parameter :: lu = 3.085678d19         ! length in kpc
    double precision, parameter :: tu = 31557600*1d6        ! time in 10^6 years
    !other constants (non-dimensionalised)
    double precision, parameter :: gcons = (6.6742d-11)*((mu*(tu**2))/(lu**3)) 
    double precision, parameter :: cms = 299792458          !light-speed in m/s
    double precision, parameter :: cs = cms*(tu/lu)
    double precision, parameter :: kap = 8d0*pi*gcons*(1d0/(cs**4))
    double precision, parameter :: kapc2 = 8d0*pi*gcons*(1d0/(cs**2))
    double precision, parameter :: Ho = (tu/(lu))*H0
    double precision, parameter :: gkr = kapc2*omega_matter*(3d0*(((Ho)**2)/(8d0*pi*gcons)))

    double precision, parameter :: lb = 3d0*omega_lambda*(((Ho)**2)/(cs*cs))
    !
    double precision, parameter :: zi = 150d0
    double precision, parameter :: zf = 0.2d0           !final redshift/position of void-lens

    ! Hamaus profile:
    double precision, parameter :: a = 2.1              !inner slope
    double precision, parameter :: b = 8.5              !outer slope
    double precision, parameter :: dc = -0.9            !central depth
    double precision, parameter :: rv = 22.1d3          !void size
    double precision, parameter :: rs = 0.88*rv         !shell size

    !initial profile
    double precision, parameter :: delta0 = -0.032      !depth
    double precision, parameter :: sigma0 =0.4          !wall thickness
    double precision, parameter :: r0 = 80              !size

    double precision, parameter :: dr = 6.25*r0/nr      !step size = prefactor*void_size/number of steps


    !decay/interaction parameters
    double precision, parameter :: vin = 0.0d0*(100/cms)    !normalised injection velocity
    double precision, parameter :: drate =1.0*(Ho/cs)       !normalised decay rate

    double precision :: ti 
    double precision :: tf 
    double precision :: dt

    !NOTE: subroutines to get and set could be contained in this module...

    !CONTAINS

    !subroutine ...

END MODULE void_parameters
