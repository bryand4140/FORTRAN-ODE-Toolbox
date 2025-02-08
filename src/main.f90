program main
    use MOD_Select_Kind, only: pv
    use MOD_ODE_Toolbox
    use MOD_ODE_Systems
    use MOD_IO_Toolbox, only: write_matrix

    implicit none

    !=======================================================
    !                    ** CONTROLS **
    logical, parameter :: export_solutions = .true.

    !=======================================================

    real(pv) :: elapsed_time
    procedure(ODE_System), pointer :: ode_ptr => null()
    real(pv) :: span(2), tolerance
    integer :: n_sys, status
    real(pv), allocatable :: IC(:), solution_matrix(:,:)
    real(pv), parameter :: pi = 3.14159265358979323846_pv


    call tic()
    print*, " " !Spacer
    print*,"==================================================="
    print*,"=============> ODE Toolbox Main <=================="

    !----------------------------------------------------------------
    ! Solve the pendulum ODE system using the RK4 Adaptive method.
    ode_ptr => Pendulum_ODE
    n_sys = 2
    span = [0.0_pv, 5.0_pv]
    tolerance = 1.0e-6_pv

    allocate(IC(n_sys))

    IC = [pi/4.0, 0.0_pv]
    
    CALL ODE_Numerical_Solve_RK4_Adaptive(span, n_sys, IC, ode_ptr,&
    solution_matrix, status, tolerance)

    if(status == 0) then
        print*,"ODE Numerical Solve RK4 Adaptive: Success"

        !Export the solutions to a csv file for plotting:
        if(export_solutions) then
            call write_matrix(solution_matrix, "ODE_Pendulum_Solution.csv",scientific=.true.)
        end if
    else
        print*,"ODE Numerical Solve RK4 Adaptive: Failure"
    end if


    !----------------------------------------------------------------
    ! Solve the damped harmonic oscillator ODE system using the RK4 Adaptive method.
    ode_ptr => Damped_Harmonic_Oscillator

    n_sys = 2
    span = [0.0_pv, 25.0_pv]
    tolerance = 1.0e-6_pv
    
    IC = [1.0_pv, 0.5_pv]

    CALL ODE_Numerical_Solve_RK4_Adaptive(span, n_sys, IC, ode_ptr,&
    solution_matrix, status, tolerance)

    if(status == 0) then
        print*,"ODE Numerical Solve RK4 Adaptive: Success"

        !Export the solutions to a csv file for plotting:
        if(export_solutions) then
            call write_matrix(solution_matrix, "ODE_Damped_Harmonic_Oscillator_Solution.csv",scientific=.true.)
        end if
    else
        print*,"ODE Numerical Solve RK4 Adaptive: Failure"
    end if


    print*,"==================================================="
    print*," " !Spacer
    elapsed_time = toc()
    if(elapsed_time > 60.0_pv) then
        write(*,"(A, ES12.4, A)") "Elapsed time = ", elapsed_time/60.0_pv, " minutes"
    else
        write(*,"(A, ES12.4, A)") "Elapsed time = ", elapsed_time, " seconds"
    end if

    contains



end program main