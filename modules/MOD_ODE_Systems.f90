module MOD_ODE_Systems
    use MOD_Select_Kind, only: pv

    implicit none

    private


contains

    ! BLasius SYS of 1st ORDER EQNS
    subroutine Blasius_Momentum_Equation(n, t, x , xdot)
        !Blasius Solution
        integer, intent(in)   :: n
        real(pv), intent(in)  :: t
        real(pv), intent(in)  :: x(n)
        real(pv), intent(out) :: xdot(n)

        xdot(1) = x(2)
        xdot(2) = x(3)
        xdot(3) = -0.50d0 * x(1) * x(3)
    end subroutine Blasius_Momentum_Equation



end module MOD_ODE_Systems