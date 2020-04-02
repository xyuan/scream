module scream_history
!----------------------------------------------------------------------
! The scream history module provides the user interface for SCREAM's
! history output capabilities.
!---------------------------------------------------------------------

  use pio,  only: pio_createfile, PIO_CLOBBER

  -------------------------------------------------------------------
  subroutine create_history_file(fname)
    character(len=*), intent(in) :: fname
    integer                      :: ierr
    integer                      :: mode

    mode = PIO_CLOBBER

    ierr = pio_createfile( , , ,fname,mode)
  end create_history_file

end module scream_history
