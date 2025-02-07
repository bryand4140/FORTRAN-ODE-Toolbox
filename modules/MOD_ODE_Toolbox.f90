module MOD_ODE_Toolbox
    use MOD_Select_Kind, only: pv

    implicit none

    integer, public :: start_time
    integer, public :: count_rate
    public :: tic, toc


contains


    ! BLasius SYS of 1st ORDER EQNS
    subroutine system_of_odes(t,x,SV)
        !Blasius Solution
        real(pv), intent(in)  :: t
        real(pv), intent(in)  :: x(3)
        real(pv), intent(out) :: SV(3)

        SV(1) = x(2)
        SV(2) = x(3)
        SV(3) = -0.50d0 * x(1) * x(3)
    end subroutine system_of_odes




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