program main
    use MOD_Select_Kind, only: pv
    use MOD_ODE_Toolbox
    use MOD_ODE_Systems

    implicit none

    real(pv) :: elapsed_time

    call tic()
    print*, " " !Spacer
    print*,"==================================================="
    print*,"=============> ODE Toolbox Main <=================="






    print*,"==================================================="
    print*," " !Spacer
    elapsed_time = toc()
    if(elapsed_time > 60.0_pv) then
        write(*,"(A, ES15.4)") "Elapsed time = ", elapsed_time/60.0_pv, " minutes"
    else
        write(*,"(A, ES15.4)") "Elapsed time = ", elapsed_time, " seconds"
    end if

    contains



end program main