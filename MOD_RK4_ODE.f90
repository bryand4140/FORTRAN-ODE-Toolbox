module MOD_RK4_ODE
    !Author: Bryan Durham
    !Module Description: This module calculates the flat plate solutions to the LD equations
    ! and uses a shooting method to find the initial conditions for a given edge Mach number,
    ! Prandtl Number, and value of cp/cv (gamma). 

    !Last Update: 12/13/2023

    implicit none

    ! Define a kind type parameter:
    !   p = 6  = real(4)
    !   p = 15 = real(8)
    !   p = 33 = real(16)
    integer, parameter :: pv = selected_real_kind(p = 15)

contains
    !LD - FLAT PLATE SYS OF 1st ORDER EQNS
    subroutine system_of_odes(t, x, gamma, Me, Pr, SV)
        !LD - FLAT PLATE SYSTEM OF ODEs
        real(pv), intent(in)  :: t
        real(pv), intent(in)  :: x(6)
        real(pv), intent(in)  :: gamma  !Ratio of specific heats
        real(pv), intent(in)  :: Me     !Mach Number at the edge of the BL
        real(pv), intent(in)  :: Pr     !Prandtl Number
        real(pv), intent(out) :: SV(6)

        !Local Variables
        real(pv) :: n, C, CP, g

        g = x(4)
        n  = 0.72;
        C  = g**(n-1);
        CP = x(5)*(n-1)*g**(n-2);
        

        SV(1) = x(2); 
        SV(2) = x(3); 
        SV(3) = -(x(3)/C)*(CP + x(1));
        SV(4) = x(5);
        SV(5) = -(Pr/C) * (CP*x(5)/Pr + x(1)*x(5) + C * (gamma - 1) * Me**2 * (x(3))**2);
        SV(6) = sqrt(2.0d0)*x(4);

    end subroutine system_of_odes

    ! ! BLasius SYS of 1st ORDER EQNS
    ! subroutine system_of_odes(t,x,SV)
    !     !Blasius Solution
    !     real(pv), intent(in)  :: t
    !     real(pv), intent(in)  :: x(3)
    !     real(pv), intent(out) :: SV(3)

    !     SV(1) = x(2)
    !     SV(2) = x(3)
    !     SV(3) = -0.50d0 * x(1) * x(3)
    ! end subroutine system_of_odes

! subroutine generate_LD_flatplate_sols(indp_var_array, n_sys, initial_conds, gamma, Me, Pr, X_SOL)
!     INTEGER, INTENT(IN)  :: n_sys !Number of equations in the system of first order ODEs
!     real(pv), intent(in)  :: indp_var_array(:) !Independant variable array
!     real(pv), intent(in)  :: initial_conds(n_sys) !Initial conditions
!     real(pv), intent(in)  :: gamma  !Ratio of specific heats
!     real(pv), intent(in)  :: Me     !Mach Number at the edge of the BL
!     real(pv), intent(in)  :: Pr     !Prandtl Number
!     real(pv), intent(out) :: X_SOL(size(indp_var_array, DIM = 1), n_sys) !Matrix to store the solutions
   

! end subroutine generate_LD_flatplate_sols


!---------------------------------------------------------------------------------------------------------
!                                     **ODE SOLVERS **

subroutine ODE_RK4(indp_var_array, n_sys, initial_conds, gamma, Me, Pr, solution)
    ! Inputs:
    !   indp_var_array - Array of independent variables (e.g., time points)
    !   initial_conds   - Initial conditions of the state variables
    !   n_sys          - The number of equations in the system of ODEs

    ! Outputs:
    !   solution  - Matrix containing the solution(s)
    implicit none

    INTEGER, INTENT(IN)  :: n_sys !Number of equations in the system of first order ODEs
    real(pv), intent(in),DIMENSION(:)  :: indp_var_array !Independant variable array
    real(pv), intent(in)  :: initial_conds(n_sys) !Initial conditions
    real(pv), intent(in)  :: gamma   !Ratio of specific heats
    real(pv), intent(in)  :: Me      !Mach Number at the edge of the BL
    real(pv), intent(in)  :: Pr      !Prandtl Number
    real(pv), intent(out),ALLOCATABLE :: solution(:,:) !Matrix to store the solutions
   
    integer :: i, n
    real(pv) :: h, t
    real(pv), ALLOCATABLE :: x(:), k1(:), k2(:), k3(:), k4(:)

    ALLOCATE(x(n_sys), k1(n_sys), k2(n_sys), k3(n_sys), k4(n_sys))
    n = size(indp_var_array, DIM = 1)  ! Number of points in the independent variable array

    ALLOCATE(solution(n,n_sys))

    ! Initialize the first row of the solution matrix with initial conditions
    solution(1, :) = initial_conds
    x = initial_conds

    ! Iteratively apply RK4
    do i = 1, n - 1
        h = indp_var_array(i+1) - indp_var_array(i)  ! Step size
        t = indp_var_array(i)

        ! RK4 calculations
        call system_of_odes(t, x, gamma, Me, Pr, k1)
        call system_of_odes(t + h/2.d0, x + h/2.d0*k1, gamma, Me, Pr, k2)
        call system_of_odes(t + h/2.d0, x + h/2.d0*k2, gamma, Me, Pr, k3)
        call system_of_odes(t + h, x + h*k3, gamma, Me, Pr, k4)

        !Calculate the updated step
        x = x + h/6.d0*(k1 + 2.d0*k2 + 2.d0*k3 + k4)

        !Store the current values in the solution:
        solution(i+1, :) = x
    end do

    deallocate(x, k1, k2, k3, k4)
end subroutine ODE_RK4


SUBROUTINE linspace(x_start, x_end, num_points, array)
    ! Arguments
    real(pv), INTENT(IN) :: x_start, x_end  ! x_start and end values
    INTEGER, INTENT(IN) :: num_points       ! Number of points
    real(pv), DIMENSION(num_points), INTENT(OUT) :: array  ! Output array

    ! Local variables
    INTEGER :: i
    real(pv) :: step

    IF (num_points > 1) THEN
        step = (x_end - x_start) / REAL(num_points - 1)
        DO i = 1, num_points
            array(i) = x_start + REAL(i - 1) * step
        END DO
    ELSE IF (num_points == 1) THEN
        array(1) = x_start
    ELSE
        PRINT *, 'Error: num_points must be >= 1'
    END IF
END SUBROUTINE linspace

end module MOD_RK4_ODE