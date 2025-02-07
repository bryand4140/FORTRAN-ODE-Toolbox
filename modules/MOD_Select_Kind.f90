module MOD_Select_Kind
    implicit none

    ! The selected_real_kind function is used to define a real variable with a specific precision.
    ! The first argument is the number of decimal digits of precision, and the second argument is the range of the exponent.
    integer, parameter, public :: pv = selected_real_kind(15, 307)


end module MOD_Select_Kind