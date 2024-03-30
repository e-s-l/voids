!
MODULE program_parameters
    implicit none

    integer, parameter :: nx = 12           !number of variables in the system of equations
    !SPACE & TIME DOMAIN RESOLUTION:
    integer, parameter :: nr = 200          !number of radial grid points
    integer, parameter :: nt = 500          !number of time grid points, might make this variable later. 
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

    double precision :: X(nr,nx)        !array of variables across space   
    double precision :: Xini(nr,nx)        !initial  
    double precision :: rad(nr)         !vector of radius values for each grid point (depends on model choice)
    
    ! for plotting purposes
    integer :: counter = 0
    
    !!!!!!!!!!!!!!!!!!!!
    !FOR PROGRAM CONTROL
    !!!!!!!!!!!!!!!!!!!!

    !SET WHICH KIND OF EVOLUTION TO BE PREFORMED:
    logical, parameter ::  lrs_evolution_mode = .false. 
    logical, parameter ::  ltb_evolution_mode = .true. 
    logical, parameter ::  linear_evolution_mode = .true. 

    !
    character(len=3) :: mode   

    !debug:
    logical, parameter ::  write_data = .true. 
    logical, parameter ::  debug_mode = .true. 

    !!!!!!!!!!!!!!!!!!!
    ! INITIALISATIONS:
    !!!!!!!!!!!!!!!!!!!



    call timelcdm(zi, ti)
    call timelcdm(zf, tf)

    dt = (tf - ti) / (1.0*nt)
    !!!
    if (debug_mode) then
        print *, 'DEBUG: dt = ', dt
    endif
    !!!


    !!!!!!!!!!!!!!!
    call initial(rad, Xini)
    !!!!!!!!!!!!!!!

  !  print *, "gkr, lb = ", gkr, lb
  !  print *, X(:,1)
    !
    if (write_data) then
        mode="ini"
        call write_delta_to_file(mode, counter, rad, Xini)      !need to redimensionalise the outputs!
    !    call write_X_to_file(counter, X)
    endif

    !!!!!!!!!!!
    !EVOLUTION
    !!!!!!!!!!!

    !
    if (lrs_evolution_mode) then
        X = Xini
        mode="LRS"
         call lrs_evolution(rad,X)  !,R_lrs) ?
         if (write_data) then
            call write_delta_to_file(mode, counter, rad, X)  
        endif

    endif
    !
    if (ltb_evolution_mode) then
        X = Xini
        mode="LTB"
        call ltb_evolution(rad,X)
        if (write_data) then
            call write_delta_to_file(mode, counter, rad, X)  
        endif
    endif
    !
    if (linear_evolution_mode) then
        X = Xini
        mode = "LIN"   
        call linear_evolution(rad, X)
        if (write_data) then
            call write_delta_to_file(mode, counter, rad, X)  
        endif
    endif 

    !
    !!!!!!!!!!!!!
    print *, ':)'
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CONTAINS
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE linear_evolution(rad, X)
            implicit none

            double precision, dimension(nr) ::  rad
            double precision, dimension(nr,nx) ::  X, Xint     
            double precision :: DD0, a_flrw, delta
            integer :: I = 1
        
            call Growth(zi,DD0)  
            a_flrw = (1.0+zi)/(1.0+0.0d0)
        !    X = Xint 

            do I= 1,nr
                delta = ((X(I,1)/(gkr*(1.0d0 + zi)**3))-1.0)/DD0
                X(I,1) = (1.0+ delta)*gkr*a_flrw**3
                rad(I) = dr*I*a_flrw*1d-3
            enddo

        END SUBROUTINE

        
        SUBROUTINE Growth(z,DD0)
            implicit none
            double precision :: z,Om,Ol,DD0,z_equality,hp,D
            double precision :: omega,omegal,oz,olz,theta_cmb

            Om = omega_matter
            Ol = omega_lambda
            hp = kapc2*omega_matter

            omega=Om
            omegal = Ol
            theta_cmb = (2.728/2.7)
            z_equality = -1.0 + (2.5*1d4*Om*hp*hp*(1.0/(theta_cmb**4)))
            oz=omega*(1.0+z)**3/(omegal+(1.0-omegal-omega)*(1.0+z)**2+omega*(1.0+z)**3)
            olz=omegal/(omegal+(1.0-omegal-omega)*(1.0+z)**2.0+omega*(1.0+z)**3)
            D = real((1.0+z_equality)/(1.0+z)*5.0*oz/2.0*(oz**(4.0/7.0)-olz+(1.0+oz/2.0)*(1.0+olz/70.0))**(-1.0))
            DD0= 1.0*D/real((1.0+z_equality)*5.0*omega/2.0*(omega**(4.0/7.0)-omegal+(1.0+omega/2.0)*(1.0+omegal/70.0))**(-1.0))

        END SUBROUTINE 

!------

        SUBROUTINE ltb_evolution(R_ltb,X)
            implicit none 

            integer :: I
            double precision, dimension(nr,nx) ::  X
            double precision :: rho,tht,shr,wey,m,mr,k,kr
            double precision :: re,rre,r,rr,rt,rtr, R_ltb(nr)

            do I=1,nr
                r = I*dr
                call ltb_initial(r,rho,tht,shr,wey,m,mr,k,kr)

                call ltb_shell_evolution(r,m,mr,k,kr,re,rre)

                r = re
                rr = rre
                rt = dsqrt(2.0d0*(m/r) - k + (lb/3.0d0)*r*r)
                rtr = ((mr/r) - (m/(r**2))*rr - 0.5d0*kr + (lb/3d0)*r*rr)/rt
            !- density
                rho = 2d0*(mr)/(r*r*(rr))
            !- expansion
                tht =  (rtr + 2.0*rt*(rr/r))/(rr )
            !- shear 
                shr = -(2.0/3.0)*dsqrt( (((rtr - rt*(rr/r))/(rr ) )**2) )
            !- weyl 
                wey =  -2*((m/r**3) - (1d0/6.0)*rho)   

                X(I,:) = [rho,0.0d0,tht,shr,wey,1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0]
            
                R_ltb(I) = re*1d-3

            enddo


        END SUBROUTINE ltb_evolution

                !----------------------------------------------------------------
        subroutine ltb_initial(r,rho,tht,shr,wey,m,mr,k,kr)
        implicit none

        logical :: compensated 

        double precision Am,Ak,r0,dlr,amp,ho,Om
        double precision dta0, dta1, dta2
        double precision d0el,d1el,d2el
        double precision m,mr,mrr
        double precision k,kr,krr,cti
        double precision r,rr,rt,rtr
        double precision rho,tht,shr,wey,epr

        !compensated = .false.
        compensated = .true.

    !    gkr=cpar(5)*cpar(18)
    !    lb = cpar(7)

        cti = ti
        epr = 0.0d0

        Om = omega_matter
        amp = delta0
        ho = Ho/cs
        dlr = sigma0*r0

     
            
        if(compensated) then 
            Am = (1.0d0/6.0d0)*gkr
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
        endif
        r = r
        rr = 1.0
        rt = dsqrt(2.0d0*(m/r) - k + (lb/3.0d0)*r*r)
        rtr = ((mr/r) - (m/(r**2))*rr - 0.5d0*kr + (lb/3d0)*r*rr)/rt
        !- density 
        rho = 2d0*(mr)/(r*r*(rr))
        !- expansion 
        tht =  (rtr + 2.0*rt*(rr/r))/(rr)
        !- shear 
        shr = -(2.0/3.0)*dsqrt( (((rtr - rt*(rr/r))/(rr) )**2) )
        !- weyl   
        wey =  -2*((m/r**3) - (1d0/6.0)*rho)


        end subroutine ltb_initial

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine timelcdm(zo, ctt) 
            implicit none

            double precision zo,rhb,rhzo,x,arsinh,tzo,ctt
            double precision ztt,thb

       !     print *, "gkr, lb = ", gkr, lb

            rhzo = gkr*((1.0d0+zo)**3 )
            x = dsqrt(lb/rhzo)
            arsinh = dlog(x + dsqrt(x*x + 1d0))
            ctt = (dsqrt((4d0)/(3d0*lb)))*arsinh

        !    print *, "z,t = ", zo, ctt

        end subroutine timelcdm

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine ltb_shell_evolution(r,m,mr,k,kr,re,rre)
            implicit none

            integer :: I
            double precision :: m,mr,k,kr
            double precision :: re,rre,ri, rri
            double precision :: r,rt,rtr,lb
            double precision :: rk1,rk2,rk3,rk4
            double precision :: rr1,rr2,rr3,rr4

            ri = r
            rri = 1.0d0
           
            do I=1,nt
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
            integer, intent(inout) :: counter
            double precision, dimension(nr,nx), intent(inout) :: X

            write(filename, '(A,:,(I1), A)') 'X.lrs.',counter,'.dat' 
            open(666, file=filename, status="replace")
            write(666,*) X

            counter = counter+1

        end subroutine

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine write_delta_to_file(mode, counter, rad, X)
            implicit none

            character(len=50) filename
            integer, intent(inout) :: counter
            double precision, dimension(nr,nx), intent(inout) :: X
            double precision, dimension(nr), intent(inout) :: rad
            integer :: I
            character(len=3) :: mode   

            open(777, file='delta.'//mode//'.dat', status="replace")
            
        
            do I=1,nr         
                write(777,*) rad(I), (X(I,1)+X(I,2))/(X(nr-10,1)+X(nr-10,2))
            enddo    



        end subroutine

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine initial(rad, X)
            implicit none

            double precision, dimension(nr,nx), intent(inout) :: X
            double precision, dimension(nr) :: rad
            double precision  :: rho,tht,shr,wey
            double precision :: r, zo, zz, tevo, ctf, cto
            double precision :: a_flrw
            integer :: I

            ! initial values: redshift, time instant, density, and expansion rate (the LCDM model assumed)
       
         !   call timelcdm(zi,cto)
         !   call timelcdm(zf,ctf)
         !   tevo = ctf-cto

            a_flrw = (1.0+zi)/(1.0+0.0d0)

            do I = 1,nr
                r = dr*I
                call ltb_profile(r,rho,tht,shr,wey)
                X(I,:) = [rho,0.0d0,tht,shr,wey,1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0]
                rad(I) = dr*I*a_flrw*1d-3
            enddo 

            print *, ":)"

        end subroutine

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine ltb_profile(r,rho,tht,shr,wey)
            implicit none

            double precision Am,Ak,dlr,sgm,amp,ho,Om,ran
            double precision dta0, dta1, dta2
            double precision d0el,d1el,d2el
            double precision m,mr,mrr,m2a,m2b,m0,m1,m2,m0r
            double precision k0,k,kr,krr,cti,ke
            double precision r,rr,rt,rtr
            double precision rho,tht,shr,wey,vol

            double precision rp,r1,r2,k1,k2

            !current implmentation for compensated profile only

            Om = omega_matter
            amp = delta0
            ho = Ho/cs
            dlr = sigma0*r0

            ! cpar(18) = gkr*((zi+1.0d0)**3)

          !  print *,'r =', r
            
            Am = (1.0d0/6.0d0)*gkr*((zi+1.0d0)**3)
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


            rr = 1.0
            rt = dsqrt(2.0d0*(m/r) - k + (lb/3.0d0)*r*r)
            rtr = ((mr/r) - (m/(r**2))*rr - 0.5d0*kr + (lb/3d0)*r*rr)/rt

            !- density
            rho = 2d0*(mr)/(r*r*(rr))


          !  print *,'dlr =', dlr
          !  print *,'dta1 =', dta1
          !  print *,'d1el =', d1el

            !- expansion
            tht =  (rtr + 2.0*rt*(rr/r))/(rr )
            !- shear
            shr = dsqrt( (1d0/3d0)*(((rtr - rt*(rr/r))/(rr))**2))
            !- weyl
            wey = - (m/r**3) + (1d0/6.0)*rho

        
        end subroutine

     
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END PROGRAM