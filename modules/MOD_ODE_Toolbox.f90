module MOD_ODE_Toolbox
    use MOD_Select_Kind, only: pv

    implicit none

    private

    !Declair Public Subroutines/Functions:
    integer, public :: start_time
    integer, public :: count_rate
    public :: tic, toc

    !Makes the interface available to other modules and programs:
    public :: ODE_System


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