!
MODULE program_parameters
    implicit none

    integer, parameter :: nx = 12       !number of variables in the system of equations
    !SPACE & TIME DOMAIN RESOLUTION:
    integer, parameter :: nr = 20  !number of radial grid points
    integer, parameter :: nt = 100 !number of time grid points, might make this variable later. 
    !Cosmological Model Parameters:
    double precision, parameter :: omega_matter = 0.3
    double precision, parameter :: omega_lambda = 1.0 - omega_matter
    double precision, parameter :: H0 = 67.810d0 !
    !
    double precision :: age
    double precision, parameter :: pi = 4d0*datan(1.0d0)
    !non-dimensioanlisations/units
    double precision, parameter :: mu = 1.989d45 !mass in 10^15 M_{\odot}
    double precision, parameter :: lu = 3.085678d19 !length in kpc
    double precision, parameter :: tu = 31557600*1d6 ! time in 10^6 years
    !other constants (non-dimensionalised)
    double precision, parameter :: gcons = 6.6742d-11*((mu*(tu**2))/(lu**3)) 
    double precision, parameter :: cms = 299792458 !light-speed in m/s
    double precision, parameter :: cs = cms*(tu/lu)
    double precision, parameter :: kap = 8d0*pi*gcons*(1d0/(cs**4))
    double precision, parameter :: kapc2 = 8d0*pi*gcons*(1d0/(cs**2))
    double precision, parameter :: Ho = (tu/(lu))*H0
    double precision, parameter :: gkr = kapc2*3d0*(((Ho)**2)/(8d0*pi*gcons))*omega_matter
    double precision, parameter :: lb = 3d0*omega_lambda*(((Ho)**2)/(cs*cs))
    !
    double precision, parameter :: zf = 150d0
    double precision, parameter :: zi = 0.2d0   !final redshift/position of void-lens

    ! Hamaus profile:
    double precision, parameter :: a = 2.1          !inner slope
    double precision, parameter :: b = 8.5          !outer slope
    double precision, parameter :: dc = -0.9        !central depth
    double precision, parameter :: rv = 22.1d3      !void size
    double precision, parameter :: rs = 0.88*rv     !shell size

    !initial profile
    double precision, parameter :: delta0 = -0.032 !depth
    double precision, parameter :: sigma0 =0.4     !wall thickness
    double precision, parameter :: r0 = 80         !size

    double precision, parameter :: dr = 6.25*r0/nr  !step size = prefactor*void_size/number of steps


    !decay/interaction parameters
    double precision, parameter :: vin = 0.0d0*(100/cms)  !normalised injection velocity
    double precision, parameter :: drate =1.0*(Ho/cs) !normalised decay rate

    double precision :: ti 
    double precision :: tf 

    !NOTE: the subroutines could be contained in this module rather than the main program as below...


END MODULE program_parameters

!
!

PROGRAM TwoFluidLRSEvolution
    use program_parameters

    implicit none !since don't want variables starting w/ i,j etc. immediately defined as integers!


    !!!!!!!!!!!!!!!!
    ! DECLARATIONS:
    !!!!!!!!!!!!!!!!

    !array of variables across space
    double precision :: X(nr,nx)
    double precision :: R_ltb(nr)
    
    ! for plotting purposes
    integer :: counter = 0
    
    !FOR PROGRAM CONTROL
    !SET WHICH KIND OF EVOLUTION TO BE PREFORMED:
    logical, parameter ::  lrs_evolution_mode = .false. 
    logical, parameter ::  ltb_evolution_mode = .true. 

    !!!!!!!!!!!!!!!!!!!
    ! INITIALISATIONS:
    !!!!!!!!!!!!!!!!!!!

    call initial(X)
    call write_X_to_file(counter,X)

    !!!!!!!!!!!
    !EVOLUTION
    !!!!!!!!!!!
    if (lrs_evolution_mode) then
         call lsr_evolution(X)
    else if (ltb_evolution_mode) then
        call ltb_evolution(X,R_ltb)
    else
    print *,':('
    endif

    call write_X_to_file(counter,X)





    !!!!!!!!!!!!!
    print *, ':)'
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CONTAINS
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE ltb_evolution(X,R_ltb)
            implicit none 

            integer :: I
            double precision, dimension(nr,nx) ::  X
            double precision :: rho,tht,shr,wey,m,mr,k,kr
            double precision :: re,rre,r,rr,rt,rtr, R_ltb(nr)

            do I=1,nr
                r = I*dr
            !   call ltb_initial(cpar,r,rho,tht,shr,wey,m,mr,k,kr)

                call ltb_shell_evolution(r,m,mr,k,kr,re,rre)

                r = re
                rr = rre
                rt = dsqrt(2.0d0*(m/r) - k + (lb/3.0d0)*r*r)
                rtr = ((mr/r) - (m/(r**2))*rr - 0.5d0*kr + (lb/3d0)*r*rr)/rt
            !- density
                rho = 2d0*(mr)/(r*r*(rr))
            !- expansion
                tht =  (rtr + 2.0*rt*(rr/r))/(rr )
            !- shear ! UPDATE to be consistent with Earl's convention
                shr = -(2.0/3.0)*dsqrt( (((rtr - rt*(rr/r))/(rr ) )**2) )
            !- weyl  ! UPDATE to be consistent with Earl's convention
                wey =  -2*((m/r**3) - (1d0/6.0)*rho)   

                X(I,:) = [rho,0.0d0,tht,shr,wey,1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0]
            
                R_ltb(I) = re*1d-3

            enddo


        END SUBROUTINE ltb_evolution

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine timelcdm(zo, ctt) 
            implicit none
            double precision zo,ct,lb,rhb,rhzo,x,arsinh,tzo,ti,ctt
            double precision ztt,thb
            rhzo = gkr*( (1.0d0+zo)**3 )
            x = dsqrt(lb/rhzo)
            arsinh = dlog(x + dsqrt(x*x + 1d0))
            ctt = (dsqrt((4d0)/(3d0*lb)))*arsinh
        end subroutine timelcdm

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine ltb_shell_evolution(r,m,mr,k,kr,re,rre)
            implicit none

            integer :: I, Ne
            double precision :: m,mr,k,kr
            double precision :: re,rre,ri, rri
            double precision :: r,rt,rtr,lb
            double precision :: rk1,rk2,rk3,rk4
            double precision :: rr1,rr2,rr3,rr4
            double precision :: dt, ti, tf


            Ne = 50


            call timelcdm(zi, ti)
            call timelcdm(zf, tf)

            dt = (tf - ti) / (1.0*Ne)
            ri = r
            rri = 1.0d0

            do I=1,Ne
                re = ri
                rt = dsqrt(2.0d0*(m/re) - k + (lb/3.0d0)*re*re)
                rk1 = dt*rt
                rre = rri
                rtr= ((mr/re) - (m/(re**2))*rre - 0.5d0*kr + (lb/3d0)*re*rre)/rt
                rr1 = dt*rtr
                re = ri + 0.5*rk1
                rt = dsqrt(2.0d0*(m/re) - k + (lb/3.0d0)*re*re)
                rk2 = dt*rt
                rre = rri+ 0.5*rr1
                rtr= ((mr/re) - (m/(re**2))*rre - 0.5d0*kr + (lb/3d0)*re*rre)/rt
                rr2 = dt*rtr
                re = ri + 0.5*rk2
                rt = dsqrt(2.0d0*(m/re) - k + (lb/3.0d0)*re*re)
                rk3 = dt*rt
                rre = rri+ 0.5*rr2
                rtr= ((mr/re) - (m/(re**2))*rre - 0.5d0*kr + (lb/3d0)*re*rre)/rt
                rr3 = dt*rtr
                re = ri + rk3
                rt = dsqrt(2.0d0*(m/re) - k + (lb/3.0d0)*re*re)
                rk4 = dt*rt
                rre = rri+ rr3
                rtr= ((mr/re) - (m/(re**2))*rre - 0.5d0*kr + (lb/3d0)*re*rre)/rt
                rr4 = dt*rtr
                re = ri + (rk1+2.0*(rk2+rk3)+rk4)/6.0d0
                ri = re
                rre = rri + (rr1+2.0*(rr2+rr3)+rr4)/6.0d0
                rri = rre
            enddo
            
        end subroutine ltb_shell_evolution

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine write_X_to_file(counter,X)
            implicit none

            character(len=50) filename
            
            double precision, dimension(40) :: cpar
            integer, intent(inout) :: counter
            double precision, dimension(nr,nx), intent(inout) :: X

            write(filename, '(A,:,(I1), A)') 'X.lrs.',counter,'.dat' 
            open(666, file=filename, status="replace")
            write(666,*) X

            counter = counter+1

        end subroutine

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine initial(X)
            implicit none

            double precision, dimension(nr,nx), intent(inout) :: X
            double precision  :: rho,tht,shr,wey
            double precision :: r, dr, zo, zz, zf, tevo, ctf, cto
            integer :: I

            ! initial values: redshift, time instant, density, and expansion rate (the LCDM model assumed)
       
            call timelcdm(zi,cto)
   
            call timelcdm(zf,ctf)
            tevo = ctf-cto
           
            ! cpar(23) = tevo
            ! cpar(24) = tevo/(1.0d0*nt)
            ! cpar(18) = cpar(5)*(zz**3)

            !
            do I = 1,nr
                r = dr*I

                call ltb_profile(r,rho,tht,shr,wey)
                X(I,:) = [rho,0.0d0,tht,shr,wey,1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0]

            enddo 

            print *, ":)"

        end subroutine

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine ltb_profile(r,rho,tht,shr,wey)
            implicit none


            double precision cpar(40)
            double precision Am,Ak,r0,dlr,sgm,amp,ho,Om,ran
            double precision dta0, dta1, dta2
            double precision d0el,d1el,d2el
            double precision m,mr,mrr,m2a,m2b,m0,m1,m2,m0r
            double precision k0,k,kr,krr,cti,ke
            double precision r,rr,rt,rtr,lb,ct
            double precision rho,tht,shr,wey,vol,gkr
            double precision rp,r1,r2,dr,k1,k2

            !current implmentation for compensated profile only


            gkr = cpar(18)
            lb = cpar(7)
            cti = cpar(21)
            Om = cpar(3)
            amp = cpar(17)
            ho = cpar(2)
            r0 = cpar(16)
            dlr = cpar(14)*r0

            
            Am = (1.0d0/6.0d0)*cpar(18)
            Ak = (10.0d0/3.0d0)*Am

            dta0 = dtanh((r-r0)/(2.0d0*dlr))
            dta1 = (1.0d0 - dta0**2)*(1.0d0/(2.0d0*dlr))
            dta2 = -dta0*dta1/dlr

            d0el = amp*0.5d0*(1.0d0 - dta0)
            d1el = -amp*0.5d0*dta1
            d2el = -amp*0.5d0*dta2

            m = Am*(1 + d0el)*r*r*r
            mr = Am*d1el*r*r*r + 3.0d0*(m/r)
            mrr = Am*d2el*(r**3)+3d0*Am*d1el*r*r+3d0*(mr/r)-3.0d0*(m/(r*r))

            k = Ak*d0el*r*r
            kr = Ak*d1el*r*r + 2.0d0*(k/r)
            krr = Ak*d2el*r*r+2.0d0*Ak*d1el*r+2.0d0*kr/r-2.0d0*(k/(r*r))


            r = r	
            rr = 1.0
            rt = dsqrt(2.0d0*(m/r) - k + (lb/3.0d0)*r*r)
            rtr = ((mr/r) - (m/(r**2))*rr - 0.5d0*kr + (lb/3d0)*r*rr)/rt

            !- density
            rho = 2d0*(mr)/(r*r*(rr))
            !- expansion
            tht =  (rtr + 2.0*rt*(rr/r))/(rr )
            !- shear
            shr = dsqrt( (1d0/3d0)*(((rtr - rt*(rr/r))/(rr))**2))
            !- weyl
            wey = - (m/r**3) + (1d0/6.0)*rho

        
        end subroutine

     
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END PROGRAM