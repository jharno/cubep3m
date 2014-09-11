
function nuPIDmap(p) result(ind)
    !
    ! Function that maps a partilce ID to an index in the array
    ! mass_p_nudm_fac which distinguishes between dark matter and neutrinos. 
    ! This is only used if -DNUPID is used.
    !

    implicit none

    integer(8), intent(in) :: p
    integer(4) :: ind

    ind = 2 - 1 / (1+p)

end function nuPIDmap

