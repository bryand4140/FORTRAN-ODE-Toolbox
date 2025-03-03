module MOD_ODE_Toolbox
    use MOD_Select_Kind, only: pv
    use ieee_arithmetic, only: IEEE_IS_NAN

    implicit none

    private

    !Declare Public Subroutines/Functions:
    integer, public :: start_time
    integer, public :: count_rate
    public :: tic, toc

    !Makes the interface available to other modules and programs:
    public :: ODE_System

    !Public Solvers:
    public :: ODE_Numerical_Solve_RK4
    public :: ODE_Numerical_Solve_RK4_Adaptive
    public :: ODE_Numerical_Solve_VSS
    public :: linspace


    !Interface for the different equation sets used with the SNLE Solver.
    abstract interface
        subroutine ODE_System(n, t, x, dx_dt)
            import pv
            integer, intent(in)   :: n        !Number of 1st order ODEs in the system
            real(pv), intent(in)  :: t        !Independent variable (time is used here as a placeholder)
            real(pv), intent(in)  :: x(n)     !State variable array
            real(pv), intent(out) :: dx_dt(n) !State variable derivatives
        end subroutine ODE_System
    end interface


contains
!--------------------------------------------------------------
!                    **ODE SOLVERS **

subroutine RK4_Fixed_Step(indp_var_array, n_sys, initial_conds, ODE_System_ptr, solution_matrix, status)
    !-----------------------------------------------------------------------
    ! General Description:
    !   A 4th-order Runge-Kutta integrator with a fixed "base" step size h.
    !   We integrate from indp_var_array(1) to indp_var_array(end),
    !   subdividing each interval [indp_var_array(i-1), indp_var_array(i)]
    !   into multiple steps if needed. The final sub-step is adjusted to land
    !   exactly on indp_var_array(i).
    !
    ! Inputs:
    !  - indp_var_array(:):  The array of independent-variable points where
    !       you want the solution. Must be monotonically increasing. 
    !  - n_sys:             Number of ODEs (size of state vector).
    !  - initial_conds(n_sys):  Initial condition for each state variable.
    !  - ODE_System_ptr:    Procedure pointer to a subroutine that computes
    !       dx_dt = f(t, x).  (Signature must match the abstract interface.)
    !
    ! Outputs:
    !  - solution_matrix(:,:):  A 2D array storing the solution.  Its first
    !       column is t-values, columns 2..(n_sys+1) are the states.
    !  - status:  Returns 0 if successful.
    !
    !-----------------------------------------------------------------------
    implicit none

    !-----------------------------------------------------------------------
    ! Inputs:
    real(pv), intent(in),  dimension(:) :: indp_var_array
    integer, intent(in)                 :: n_sys
    real(pv), intent(in),  dimension(n_sys) :: initial_conds
    procedure(ODE_System), pointer       :: ODE_System_ptr

    ! Outputs:
    real(pv), allocatable, intent(out) :: solution_matrix(:,:)
    integer,               intent(out) :: status

    ! Local variables
    integer :: num_steps, max_steps
    integer :: i, step
    real(pv) :: h, t, t_end, h_local
    real(pv), allocatable :: x(:)
    real(pv), allocatable :: k1(:), k2(:), k3(:), k4(:), dx_dt(:)
    real(pv), allocatable :: temp_solution(:,:)

    !-----------------------------------------------------------------------
    ! Error checking
    if (size(indp_var_array) < 2) then
        status = 2
        write(*,*) "ERROR: indp_var_array must have at least two points."
        return
    end if

    !Check to see if the elements of indp_var_array are monotonically increasing
    do i = 2, size(indp_var_array)
        if (indp_var_array(i) <= indp_var_array(i-1)) then
            status = 3
            write(*,*) "ERROR: indp_var_array must be monotonically increasing."
            return
        end if
    end do

    !-----------------------------------------------------------------------
    ! 1. Basic setup and array allocations
    !-----------------------------------------------------------------------
    num_steps = size(indp_var_array)

    ! Base step size h is determined from first two points (Assumes uniform spacing):
    h = indp_var_array(2) - indp_var_array(1)

    ! Maximum number of sub-steps = 100 per interval
    max_steps = 1 + 100*(num_steps - 1)

    allocate(x(n_sys))
    allocate(k1(n_sys), k2(n_sys), k3(n_sys), k4(n_sys), dx_dt(n_sys))

    ! Accumulate the solution in temp_solution, then trim at the end.
    allocate(temp_solution(max_steps, n_sys + 1))

    !-----------------------------------------------------------------------
    ! 2. Initialize
    !-----------------------------------------------------------------------
    x = initial_conds          ! State vector at the initial point
    step = 1                   ! First entry in row 1
    temp_solution(step, 1) = indp_var_array(1)   ! t
    temp_solution(step, 2:)   = x                ! states

    !-----------------------------------------------------------------------
    ! 3. Main loop over intervals in indp_var_array
    !-----------------------------------------------------------------------
    do i = 2, num_steps
        t     = indp_var_array(i-1)
        t_end = indp_var_array(i)

        ! Integrate from t to t_end using multiple steps of size h
        do while (t < t_end)

            ! If we're about to overshoot t_end, reduce step size
            if (t + h > t_end) then
                h_local = t_end - t
            else
                h_local = h
            end if

            ! Evaluate derivative at current state
            call ODE_System_ptr(n_sys, t, x, dx_dt)
            k1 = h_local * dx_dt

            ! Evaluate derivative at midpoint #1
            call ODE_System_ptr(n_sys, t + h_local/2._pv, x + k1/2._pv, dx_dt)
            k2 = h_local * dx_dt

            ! Evaluate derivative at midpoint #2
            call ODE_System_ptr(n_sys, t + h_local/2._pv, x + k2/2._pv, dx_dt)
            k3 = h_local * dx_dt

            ! Evaluate derivative at full step
            call ODE_System_ptr(n_sys, t + h_local, x + k3, dx_dt)
            k4 = h_local * dx_dt

            ! Update the state
            x = x + (k1 + 2._pv*k2 + 2._pv*k3 + k4) / 6._pv

            ! Advance time
            t = t + h_local

            ! Increment solution index
            step = step + 1
            if (step > max_steps) then
                ! We ran out of allocated rows => error or reallocate
                status = 1
                write(*,*) "ERROR: Exceeded maximum sub-steps = ", max_steps
                return
            end if

            ! Store results
            temp_solution(step, 1)  = t
            temp_solution(step, 2:) = x

        end do  ! while (t < t_end)
    end do  ! i=2..num_steps

    !-----------------------------------------------------------------------
    ! 4. Trim and output final solution
    !-----------------------------------------------------------------------
    ! 'step' is now the total number of stored rows
    allocate(solution_matrix(step, n_sys + 1))
    solution_matrix(:, :) = temp_solution(1:step, :)

    ! We can safely deallocate local arrays now
    deallocate(x, k1, k2, k3, k4, dx_dt, temp_solution)

    ! Indicate success
    status = 0

end subroutine RK4_Fixed_Step


SUBROUTINE RK4_Fixed_Step2(y0, y_end, h, initial_conds, ODE_System_ptr, solution_matrix, status)
    IMPLICIT NONE
    
    ! Input parameters
    REAL(pv), INTENT(IN) :: y0                    ! Initial independent variable
    REAL(pv), INTENT(IN) :: y_end                 ! Final independent variable
    REAL(pv), INTENT(IN) :: h                     ! Step size
    REAL(pv), INTENT(IN) :: initial_conds(:)      ! Initial conditions array
    PROCEDURE(ODE_System), INTENT(IN), POINTER :: ODE_System_ptr ! Pointer to the ODE system function
    
    ! Outputs
    REAL(pv), ALLOCATABLE, INTENT(OUT) :: solution_matrix(:,:)  ! Solution matrix (time steps x (system_size + 1))
    INTEGER, INTENT(OUT) :: status                ! Status flag for error handling
    
    ! Local variables
    INTEGER :: i, system_size
    REAL(pv), ALLOCATABLE :: time_points(:)
    REAL(pv), ALLOCATABLE :: y(:), k1(:), k2(:), k3(:), k4(:), temp(:)
    
    ! Initialize
    system_size = SIZE(initial_conds)
    status = 0
    
    ! Generate time points using linspace_h
    CALL linspace_h(y0, y_end, h, time_points)
    IF (.NOT. ALLOCATED(time_points)) THEN
        status = 2  ! Error in generating time points
        RETURN
    END IF
    
    ! Allocate the solution matrix first (add 1 column for independent variable)
    ALLOCATE(solution_matrix(SIZE(time_points), system_size + 1), STAT=status)
    IF (status /= 0) THEN
        status = 4  ! Memory allocation error for solution matrix
        DEALLOCATE(time_points)
        RETURN
    END IF
    
    ! Allocate working arrays
    ALLOCATE(y(system_size), k1(system_size), k2(system_size), &
             k3(system_size), k4(system_size), temp(system_size), STAT=status)
    IF (status /= 0) THEN
        status = 1  ! Memory allocation error for working arrays
        DEALLOCATE(time_points, solution_matrix)
        RETURN
    END IF
    
    ! Set initial conditions
    y = initial_conds
    solution_matrix(1,1) = time_points(1)          ! Store initial time
    solution_matrix(1,2:) = y                      ! Store initial state
    
    ! Main RK4 integration loop
    DO i = 2, SIZE(time_points)
        ! Store time point
        solution_matrix(i,1) = time_points(i)
        
        ! k1 = f(t, y)
        CALL ODE_System_ptr(system_size, time_points(i-1), y, k1)
        
        ! k2 = f(t + h/2, y + (h/2)*k1)
        temp = y + (h/2.0_pv)*k1
        CALL ODE_System_ptr(system_size, time_points(i-1) + h/2.0_pv, temp, k2)
        
        ! k3 = f(t + h/2, y + (h/2)*k2)
        temp = y + (h/2.0_pv)*k2
        CALL ODE_System_ptr(system_size, time_points(i-1) + h/2.0_pv, temp, k3)
        
        ! k4 = f(t + h, y + h*k3)
        temp = y + h*k3
        CALL ODE_System_ptr(system_size, time_points(i-1) + h, temp, k4)
        
        ! Update solution: y(n+1) = y(n) + h*(k1 + 2*k2 + 2*k3 + k4)/6
        y = y + h*(k1 + 2.0_pv*k2 + 2.0_pv*k3 + k4)/6.0_pv
        
        ! Store solution (skip first column as it's for time)
        solution_matrix(i,2:) = y
    END DO
    
    ! Clean up
    DEALLOCATE(y, k1, k2, k3, k4, temp, time_points)
    
END SUBROUTINE RK4_Fixed_Step2


subroutine ODE_Numerical_Solve_RK4(indp_var_span, n_sys, initial_conds, ODE_System_ptr, solution_matrix, status)
    !-----------------------------------------------------------------------
    ! This subroutine uses a fixed-step classical Runge-Kutta (RK4) method to
    ! integrate a system of first-order ODEs.
    !
    ! Arguments:
    !  - indp_var_span: a 2-element array giving the start and end values of
    !                   the independent variable (e.g., time = [t_start, t_end]).
    !  - n_sys        : the number of ODEs in the system.
    !  - initial_conds: an array of initial conditions (dimension n_sys).
    !  - ODE_System_ptr: a procedure pointer to a subroutine that computes the ODE
    !                   right-hand side. It must conform to the abstract interface.
    !  - solution_matrix:
    !         an allocated real array of dimension (n_sys+1, n_steps) where:
    !           * Row 1 will hold the independent variable values.
    !           * Rows 2 through n_sys+1 hold the corresponding state values.
    !  - status       : an integer output (0 for success, nonzero for failure).
    !-----------------------------------------------------------------------
    implicit none

    ! Input arguments
    real(pv), intent(in)    :: indp_var_span(2)
    integer, intent(in)     :: n_sys
    real(pv), intent(in)    :: initial_conds(n_sys)
    procedure(ODE_System), pointer, intent(in) :: ODE_System_ptr

    ! In/out argument: solution matrix should be pre-allocated by the caller.
    real(pv), intent(inout) :: solution_matrix(:, :)

    ! Output argument
    integer, intent(out)    :: status

    ! Local variables
    integer :: i, n_steps
    real(pv) :: t, dt
    real(pv), dimension(n_sys) :: x, k1, k2, k3, k4, x_temp

    ! Check that the first dimension of solution_matrix is n_sys+1.
    if ( size(solution_matrix, 1) /= n_sys + 1 ) then
        status = 1
        return
    end if

    ! Determine the number of solution points and compute the step size.
    n_steps = size(solution_matrix, 2)
    dt = ( indp_var_span(2) - indp_var_span(1) ) / real(n_steps - 1, kind=pv)

    ! Set the initial conditions into the solution matrix.
    solution_matrix(1,1)     = indp_var_span(1)
    solution_matrix(2:n_sys+1,1) = initial_conds

    ! Loop over time steps using the RK4 method.
    do i = 1, n_steps - 1
        t = solution_matrix(1,i)
        x = solution_matrix(2:n_sys+1, i)
        
        ! Compute k1 = f(t, x)
        call ODE_System_ptr(n_sys, t, x, k1)
        
        ! Compute k2 = f(t + dt/2, x + (dt/2)*k1)
        x_temp = x + 0.5_pv * dt * k1
        call ODE_System_ptr(n_sys, t + 0.5_pv * dt, x_temp, k2)
        
        ! Compute k3 = f(t + dt/2, x + (dt/2)*k2)
        x_temp = x + 0.5_pv * dt * k2
        call ODE_System_ptr(n_sys, t + 0.5_pv * dt, x_temp, k3)
        
        ! Compute k4 = f(t + dt, x + dt*k3)
        x_temp = x + dt * k3
        call ODE_System_ptr(n_sys, t + dt, x_temp, k4)
        
        ! Update the state vector:
        x = x + dt / 6.0_pv * (k1 + 2.0_pv*k2 + 2.0_pv*k3 + k4)
        
        ! Save the new time and state into the solution matrix.
        solution_matrix(1, i+1)         = t + dt
        solution_matrix(2:n_sys+1, i+1) = x
    end do

    ! Set status to success.
    status = 0
end subroutine ODE_Numerical_Solve_RK4


subroutine ODE_Numerical_Solve_RK4_Adaptive(indp_var_span, n_sys, initial_conds, ODE_System_ptr, &
    solution_matrix, status, tolerance)
    !------------------------------------------------------------
    ! Adaptive RK4 solver subroutine.
    !
    ! Arguments:
    !   indp_var_span   - Two-element array: [start, end] of independent variable.
    !   n_sys           - Number of ODEs in the system.
    !   initial_conds   - Array of initial conditions (dimension n_sys).
    !   ODE_System_ptr  - Procedure pointer to the ODE system.
    !   tolerance       - Local error tolerance for each step.
    !
    ! Outputs:
    !   solution_matrix - Allocatable solution array of size (n_steps, n_sys+1)
    !                     where column 1 holds the independent variable and
    !                     columns 2..n_sys+1 hold the state variables.
    !   status          - 0 on success
    !                   - 1 if the solution matrix is too small
    !                   - 2 if dt is non-positive
    !                   - 3 if dt becomes too small
    !------------------------------------------------------------
    implicit none
    ! Input arguments
    real(pv), intent(in)    :: indp_var_span(2)
    integer, intent(in)     :: n_sys
    real(pv), intent(in)    :: initial_conds(n_sys)
    procedure(ODE_System), pointer, intent(in) :: ODE_System_ptr
    real(pv), intent(in)    :: tolerance

    ! Output arguments
    integer, intent(out)    :: status
    real(pv), allocatable, intent(out) :: solution_matrix(:, :)

    ! Local variables
    integer :: capacity, current_step, new_capacity
    real(pv) :: t, tf, dt, dt_new, err, safety, dt_min
    real(pv), allocatable :: x(:), x_big(:), x_half(:), x_temp(:)

    ! Safety factor and a minimum allowed dt
    safety = 0.9_pv
    dt_min = 1.0e-12_pv

    ! Allocate temporary arrays of size n_sys.
    allocate(x(n_sys), x_big(n_sys), x_half(n_sys), x_temp(n_sys))

    ! Set initial conditions.
    t  = indp_var_span(1)
    tf = indp_var_span(2)
    x  = initial_conds

    ! Initial guess for dt: for example, 1/100th of the interval.
    dt = (tf - t) / 100.0_pv
    if (dt <= 0.0_pv) then
        status = 2
        return
    end if

    ! Allocate an initial solution matrix.
    capacity = 500
    ! Now the matrix has 'capacity' rows and (n_sys+1) columns.
    allocate(solution_matrix(capacity, n_sys+1))
    current_step = 1
    solution_matrix(current_step, 1)     = t
    solution_matrix(current_step, 2:n_sys+1) = x

   ! Main adaptive time stepping loop.
   do while (t < tf)
       ! Adjust dt if the next step overshoots the final time.
       if (t + dt > tf) dt = tf - t

       !--- Compute one full RK4 step over dt ---
       call RK4_Step(dt, n_sys, t, x, ODE_System_ptr, x_big)
       
       !--- Compute two half-steps (dt/2 each) ---
       call RK4_Step(dt/2.0_pv, n_sys, t, x, ODE_System_ptr, x_temp)
       call RK4_Step(dt/2.0_pv, n_sys, t + dt/2.0_pv, x_temp, ODE_System_ptr, x_half)
       
       !--- Estimate the local error ---
       err = maxval( abs( x_big - x_half ) )
       
       if (err <= tolerance) then
           ! Accept the step.
           t = t + dt
           x = x_half  ! Use the more accurate (two half-step) result.
           
           ! Increment the solution counter and save the new step.
           current_step = current_step + 1
           if (current_step > capacity) then
               ! Increase the capacity (for example, double it).
               new_capacity = capacity * 2
               call extend_solution_matrix(solution_matrix, capacity, new_capacity)
               capacity = new_capacity
           end if
           if (current_step <= size(solution_matrix, 1)) then
                solution_matrix(current_step, 1)     = t
                solution_matrix(current_step, 2:n_sys+1) = x
            else
                print *, "Error: Array bounds exceeded when storing initial conditions"
                status = 1
                return
            end if
           
           ! Compute new step size. If error is zero, try increasing dt.
           if (err == 0.0_pv) then
               dt_new = dt * 2.0_pv
           else
               dt_new = dt * safety * (tolerance / err)**(0.2_pv)
           end if
           dt = max(dt_new, dt_min)
       else
           ! Reject the step; reduce dt and try again.
           dt = dt * safety * (tolerance / err)**(0.2_pv)
           if (dt < dt_min) then
               status = 3  ! dt has become too small.
               return
           end if
       end if
   end do
   
   if (current_step < capacity) then
        block 
            real(pv), allocatable :: temp(:,:)
            allocate(temp(current_step, n_sys+1))
            temp = solution_matrix(1:current_step, :)
            call move_alloc(temp, solution_matrix)
        end block
    end if
   
   status = 0  ! Indicate success.

   contains
   subroutine RK4_Step(time_step, n_system, time, SV, ode_ptr, x_new)
        !------------------------------------------------------------
        ! Internal subroutine to take one RK4 step of size time_step.
        ! Given (time,SV) it computes x_new at time time+time_step.
        !------------------------------------------------------------
        real(pv), intent(in)    :: time_step, time
        integer, intent(in)     :: n_system
        real(pv), intent(in)    :: SV(n_system)
        procedure(ODE_System), pointer, intent(in) :: ode_ptr
        real(pv), intent(out)   :: x_new(n_system)
        
        real(pv) :: k1(n_system), k2(n_system), k3(n_system), k4(n_system)
        real(pv) :: x_temporary(n_system)
        
        call ode_ptr(n_system, time, SV, k1)
        
        x_temporary = SV + 0.5_pv * time_step * k1
        call ode_ptr(n_system, time + 0.5_pv * time_step, x_temporary, k2)
        
        x_temporary = SV + 0.5_pv * time_step * k2
        call ode_ptr(n_system, time + 0.5_pv * time_step, x_temporary, k3)
        
        x_temporary = SV + time_step * k3
        call ode_ptr(n_system, time + time_step, x_temporary, k4)
        
        x_new = SV + time_step / 6.0_pv * ( k1 + 2.0_pv * k2 + 2.0_pv * k3 + k4 )
    end subroutine RK4_Step


    subroutine extend_solution_matrix(sol_mat, old_capacity, updated_capacity)
        implicit none
        real(pv), allocatable, intent(inout) :: sol_mat(:, :)
        integer, intent(in) :: old_capacity, updated_capacity
        real(pv), allocatable :: temp(:, :)
        integer :: ncols
        
        ncols = size(sol_mat, 2)  ! Get number of columns
        allocate(temp(updated_capacity, ncols))  ! Allocate with new capacity
        temp(1:old_capacity, :) = sol_mat(1:old_capacity, :)  ! Copy existing data
        deallocate(sol_mat)
        call move_alloc(temp, sol_mat)
    end subroutine extend_solution_matrix

end subroutine ODE_Numerical_Solve_RK4_Adaptive


subroutine ODE_Numerical_Solve_VSS(indp_var_span, n_sys, initial_conds, &
    ODE_System_ptr, solution_matrix, status, tolerance, ISS, MSS)
    ! General Description:
    ! ODE Numberical Solve Variable Step Size (VSS) subroutine
    ! This subroutine solves a system of ODEs using an adaptive step size control method.
    ! It is designed to handle moderately stiff ODEs and uses the Dormand-Prince RK5 adaptive step method. 
    
    !The solver requires the following inputs:
    ! 1. indp_var_span = [x_lower, x_upper] - this array of two values specify the lower and 
    !    upper bounds for the independant variable array. The solver will determine the step 
    !    size automaticly to maintain both speed and accuracy. The independant variable array
    !    is returned as the first column in the solution vector.

    ! 2. The initial conditions (initial_conds) are the values of the state variables at the
    !    starting point of the solution. The number of initial conditions must match the number
    !    of equations in the first order system (size(IC) = n_sys).

    ! 3. The number of equations in the system (n_sys) defines the size of the state vector and
    !    the number of equations in the system of ODEs. The size of the output from the ODE_System
    !    subroutine must match the number of equations in the system.

    ! 4. The solution_matrix (solution) stores the values of the state variables at each step of
    !    the solution. The size of the solution matrix is determined dynamically. It should be
    !    initialized as

        ! >>> real(pv), allocatable :: solution(:,:)
    
    ! 5. If the unknown variables are x1, x2, ..., xn, then the solution matrix will contain
    !    the values of  x1, x2, ..., xn at each step in the folling matrix form.

    !                         [t(0),  x1(0), ,x2(0), ..., xn(0)]
    !    Solution Matrix ==   [t(1),  x1(1), ,x2(1), ..., xn(1)]
    !                         [t(2),  x1(2), ,x2(2), ..., xn(2)]
    !
    !    Note that the solution matrix does NOT contain the independant variable array.

    ! 6. The status variable (status) indicates the status of the solution. If the solver runs
    !    correctly without errors, the status is set to 0. If any NaN values are encountered in
    !    the solution, the status is set to 1 to alert the user of the issue. For an invalid 
    !    solver option, the status is set to 2.

    ! 7. The ODE_System_ptr is a procedure pointer to the ODE_System subroutine that defines the
    !    system of ODEs to be solved. The ODE_System subroutine must be implemented by the user
    !    and should define the system of ODEs to be solved.

    !NOTE on ODE_System subroutine:
    !The ODE_System subroutine must be implemented by the user and should define the system of
    !ODEs to be solved. The subroutine takes the independent variable, the state vector (x),
    !and the derivative vector (dx_deta) as inputs. The derivative vector dx_deta should contain
    !the derivatives of the state variables with respect to the independent variable. The size of
    !the state vector x and the derivative vector dx_deta must match the number of equations in
    !the system (n_sys). 

    implicit none

    ! Input arguments
    real(pv), intent(in), dimension(2) :: indp_var_span  ! Independent variable span [x0, x_end]
    integer, intent(in) :: n_sys                          ! Number of equations in the system
    real(pv), intent(in), dimension(n_sys) :: initial_conds  ! Initial conditions
    procedure(ODE_System), pointer :: ODE_System_ptr      ! Pointer to the ODE system subroutine

    !Optional arguments
    real(pv), intent(in), optional :: tolerance   ! Error tolerance
    real(pv), intent(in), optional :: ISS         !Initial step size
    real(pv), intent(in), optional :: MSS         !Maximum step size

    ! Output arguments
    real(pv), allocatable, intent(out) :: solution_matrix(:,:)  ! Solution matrix [t, y1, y2, ..., yn]
    integer, intent(out) :: status                             ! Status of the solution

    ! Local variables
    real(pv), parameter :: h_min = 1.0e-10_pv      ! Minimum step size
    real(pv), parameter :: safety_factor = 0.9_pv  ! Safety factor for step size control
    integer, parameter  :: max_step = 100000       ! Maximum number of steps

    real(pv) :: tol
    real(pv) :: t, t_end, h, h_new, err, h_max
    real(pv), allocatable :: x(:), y_high(:), y_low(:)
    real(pv), allocatable :: k1(:), k2(:), k3(:), k4(:), k5(:), k6(:), k7(:)
    real(pv), allocatable :: indp_var_array(:)
    real(pv), allocatable :: solution(:,:)
    integer :: step

    !Set the maximum step size if not present
    if (present(MSS)) then
        h_max = MSS
    else
        h_max = 1.0e-2_pv
    end if

    ! Set the tolerance if not present:
    if (present(tolerance)) then
        tol = tolerance
    else
        tol = 1.0e-6_pv
    end if

    ! Set the initial step size if not present:
    if (present(ISS)) then
        h = ISS
    else
        h = (indp_var_span(2) - indp_var_span(1)) / 1000.0_pv
    end if

    ! Initialize variables for the start and end of the interval
    t = indp_var_span(1)
    t_end = indp_var_span(2)

    ! Allocate memory for the state vector and the k values
    allocate(x(n_sys))
    allocate(y_high(n_sys), y_low(n_sys))
    allocate(k1(n_sys), k2(n_sys), k3(n_sys), k4(n_sys), k5(n_sys), k6(n_sys), k7(n_sys))
    allocate(solution(1000, n_sys))
    allocate(indp_var_array(1000))

    x = initial_conds
    solution(1, :) = x
    indp_var_array(1) = t

    status = 0  ! Initialize status
    step = 1

    do while (t < t_end .and. step < max_step)


        ! Ensure h does not exceed the remaining interval
        if (t + h > t_end) then
            h = t_end - t
        end if

        ! Compute k1 to k7
        call ODE_System_ptr(n_sys, t, x, k1)

        call ODE_System_ptr(n_sys, t + 1.0_pv/5.0_pv*h, x + h*(1.0_pv/5.0_pv)*k1, k2)

        call ODE_System_ptr(n_sys, t + 3.0_pv/10.0_pv*h, x + h*(3.0_pv/40.0_pv*k1 + 9.0_pv/40.0_pv*k2), k3)

        call ODE_System_ptr(n_sys, t + 4.0_pv/5.0_pv*h, x + h*(44.0_pv/45.0_pv*k1 - 56.0_pv/15.0_pv*k2 + 32.0_pv/9.0_pv*k3), k4)

        call ODE_System_ptr(n_sys, t + 8.0_pv/9.0_pv*h, x + h*(19372.0_pv/6561.0_pv*k1 - 25360.0_pv/2187.0_pv*k2 + &
                         64448.0_pv/6561.0_pv*k3 - 212.0_pv/729.0_pv*k4), k5)

        call ODE_System_ptr(n_sys, t + h, x + h*(9017.0_pv/3168.0_pv*k1 - 355.0_pv/33.0_pv*k2 + &
                         46732.0_pv/5247.0_pv*k3 + 49.0_pv/176.0_pv*k4 - 5103.0_pv/18656.0_pv*k5), k6)

        call ODE_System_ptr(n_sys, t + h, x + h*(35.0_pv/384.0_pv*k1 + 500.0_pv/1113.0_pv*k3 + &
                         125.0_pv/192.0_pv*k4 - 2187.0_pv/6784.0_pv*k5 + 11.0_pv/84.0_pv*k6), k7)

        ! High-order solution (5th order)
        y_high = x + h * (35.0_pv/384.0_pv*k1 + 0.0_pv*k2 + 500.0_pv/1113.0_pv*k3 + &
                          125.0_pv/192.0_pv*k4 - 2187.0_pv/6784.0_pv*k5 + 11.0_pv/84.0_pv*k6)

        ! Low-order solution (4th order)
        y_low = x + h * (5179.0_pv/57600.0_pv*k1 + 0.0_pv*k2 + 7571.0_pv/16695.0_pv*k3 + &
                         393.0_pv/640.0_pv*k4 - 92097.0_pv/339200.0_pv*k5 + 187.0_pv/2100.0_pv*k6 + 1.0_pv/40.0_pv*k7)

        ! Compute the error estimate
        err = maxval(abs(y_high - y_low))

        ! Adaptive step size control
        if (err <= tol) then
            ! Accept the step
            t = t + h
            x = y_high
            step = step + 1

            ! Store the results
            if (step > size(solution, 1)) then
                ! Expand the arrays
                call expand_arrays(solution, indp_var_array)
            end if

            indp_var_array(step) = t
            solution(step, :) = x
        else
            ! Reject the step
            ! Do not update t or x
        end if

        ! Compute new step size
        if (err == 0.0_pv) then
            ! Prevent division by zero
            h_new = h_max
        else
            h_new = safety_factor * h * (tol / err)**(1.0_pv / 5.0_pv)
            h_new = min(max(h_new, h_min), h_max)
        end if

        h = h_new

        ! Check for NaN values
        if (any(IEEE_IS_NAN(x))) then
            print *, 'Error: NaN encountered in computations at time t = ', t
            status = 1  ! Set status to 1 to alert the user of the issue
            exit
        end if

    end do

    if (step >= max_step) then
        print *, 'Warning: Maximum number of steps reached.'
        status = 2  ! Set status to 2 to indicate maximum steps reached
    end if

    ! Trim arrays to actual size
    allocate(solution_matrix(step, n_sys + 1))
    solution_matrix(:, 1)  = indp_var_array(1:step)
    solution_matrix(:, 2:) = solution(1:step, :)

    ! Deallocate arrays
    deallocate(x, y_high, y_low, k1, k2, k3, k4, k5, k6, k7)
    deallocate(solution, indp_var_array)

end subroutine ODE_Numerical_Solve_VSS

!--------------------------------------------------------------
!                 ** Helper Subroutines **

subroutine expand_arrays(solution, indp_var_array)
    ! General Description: This subroutine expands the size of the solution matrix and the independent variable array.
    ! The new size is twice the current size.
    ! The expanded elements are initialized to zero.
    ! Example: If the original size is 1000, the new size will be 2000.

    implicit none
    real(pv), allocatable, intent(inout) :: solution(:,:), indp_var_array(:)
    integer :: new_size

    new_size = size(solution, 1) * 2

    call resize_array(solution, new_size)
    call resize_vector(indp_var_array, new_size)
end subroutine expand_arrays


subroutine resize_array(array, new_size)
    ! General Description: This subroutine resizes a 2D array to a new size.
    ! Example: If the original array is [[1, 2], [3, 4]] and the new size is (3, 2),
    ! the resized array will be [[1, 2], [3, 4], [0, 0]].
    
    implicit none
    real(pv), allocatable, intent(inout) :: array(:,:)
    integer, intent(in) :: new_size
    real(pv), allocatable :: temp(:,:)
    integer :: old_size

    old_size = size(array, 1)
    allocate(temp(new_size, size(array, 2)))
    temp(1:old_size, :) = array
    deallocate(array)
    array = temp
end subroutine resize_array


subroutine resize_vector(vector, new_size)
    ! General Description: This subroutine resizes a 1D vector to a new size.
    ! Example: If the original vector is [1, 2, 3] and the new size is 5,
    ! the resized vector will be [1, 2, 3, 0, 0].

    implicit none
    real(pv), allocatable, intent(inout) :: vector(:)
    integer, intent(in) :: new_size
    real(pv), allocatable :: temp(:)
    integer :: old_size

    old_size = size(vector)
    allocate(temp(new_size))
    temp(1:old_size) = vector
    deallocate(vector)
    vector = temp
end subroutine resize_vector

!--------------------------------------------------------------
!          ** Useful Array Generation Subroutines **

SUBROUTINE linspace(x_start, x_end, array)
    ! Generates a linearly spaced array of real numbers between x_start and x_end
    ! with the number of elements inferred from the size of the array.
    ! Inputs:
    !  x_start = Starting value
    !  x_end = Ending value
    ! Outputs:
    !  array = Linearly spaced array of real numbers

    implicit none
    real(pv), INTENT(IN) :: x_start, x_end  ! x_start and end values
    real(pv), DIMENSION(:), INTENT(OUT) :: array  ! Output array

    ! Local variables
    INTEGER :: i, num_points
    real(pv) :: step

    num_points = SIZE(array)
    
    IF (num_points > 1) THEN
        step = (x_end - x_start) / REAL(num_points - 1)
        DO i = 1, num_points
            array(i) = x_start + REAL(i - 1) * step
        END DO
    ELSE IF (num_points == 1) THEN
        array(1) = x_start
    ELSE
        PRINT *, 'Error: array size must be >= 1'
    END IF
END SUBROUTINE linspace


!---------------------------------------------------------------
!                 ** Timing Subroutines **

subroutine tic()
    call system_clock(count_rate=count_rate)
    call system_clock(start_time)
end subroutine tic


function toc() result(elapsed_time)
    real(pv) :: elapsed_time
    integer :: finish_time

    call system_clock(finish_time)
    elapsed_time = real(finish_time - start_time, pv) / real(count_rate, pv)
end function toc


end module MOD_ODE_Toolbox
