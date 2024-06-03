MODULE void_subroutines  

    use void_parameters
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

                !!!
         !   if (debug_mode) then
         !       print *, 'DEBUG: dt = ', dt
         !        print *, 'DEBUG: nt = ', nt
         !   endif
            !!!


           
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

END MODULE void_subroutines  
