!
!
!
!

PROGRAM voids

    ! include the module files (parameters, subroutines)

    use void_parameters
    use void_subroutines

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
    logical, parameter ::  ltb_evolution_mode = .false. 
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

    !!!!!!!!!!!!!
    print *, 'DEBUG: Main program complete :)'
    !!!!!!!!!!!!!

END PROGRAM voids
