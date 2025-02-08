module MOD_ODE_Systems
    use MOD_Select_Kind, only: pv

    implicit none

    private

    public :: Blasius_Momentum_Equation
    public :: Pendulum_ODE
    public :: Damped_Harmonic_Oscillator


contains


subroutine Blasius_Momentum_Equation(n, t, x , xdot)
    !Blasius equations for a flat place CPG boundary layer.
    integer, intent(in)   :: n
    real(pv), intent(in)  :: t
    real(pv), intent(in)  :: x(n)
    real(pv), intent(out) :: xdot(n)

    xdot(1) = x(2)
    xdot(2) = x(3)
    xdot(3) = -0.50d0 * x(1) * x(3)
end subroutine Blasius_Momentum_Equation


subroutine Pendulum_ODE(n, t, x, xdot)
    ! ODE system for a simple pendulum (small or large amplitude).
    ! The pendulum equation: θ'' + (g/L)*sin(θ) = 0 is transformed to a first-order system:
    ! x(1) = θ     and     x(2) = dθ/dt.
    ! Thus, we have:
    !   dθ/dt    = x(2)
    !   d²θ/dt² = - (g/L) * sin(x(1))
    integer, intent(in)   :: n
    real(pv), intent(in)  :: t
    real(pv), intent(in)  :: x(n)
    real(pv), intent(out) :: xdot(n)

    ! Define constants for gravity and pendulum length.
    real(pv), parameter :: g = 9.81_pv
    real(pv), parameter :: L = 1.0_pv

    if (n < 2) then
        return  ! Ensure that at least two state variables are provided.
    end if

    xdot(1) = x(2)
    xdot(2) = - (g/L) * sin(x(1))
end subroutine Pendulum_ODE


subroutine Damped_Harmonic_Oscillator(n, t, x, xdot)
    ! ODE system for a damped harmonic oscillator.
    ! The equation: x'' + 2ζω x' + ω² x = 0 is transformed to a first-order system:
    ! x(1) = x     and     x(2) = dx/dt.
    ! Thus:
    !   dx/dt    = x(2)
    !   d²x/dt² = -2ζω x(2) - ω² x(1)
    integer, intent(in)   :: n
    real(pv), intent(in)  :: t
    real(pv), intent(in)  :: x(n)
    real(pv), intent(out) :: xdot(n)

    ! Define constants for the oscillator:
    real(pv), parameter :: zeta  = 0.3_pv    ! Damping ratio
    real(pv), parameter :: omega = 1.0_pv    ! Natural frequency

    if (n < 2) then
        print*, "Error: At least two state variables are required."
        return  ! Ensure that at least two state variables are provided.
    end if

    xdot(1) = x(2)
    xdot(2) = - 2.0_pv * zeta * omega * x(2) - omega**2 * x(1)
end subroutine Damped_Harmonic_Oscillator



end module MOD_ODE_Systems