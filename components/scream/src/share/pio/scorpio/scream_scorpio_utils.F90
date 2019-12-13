#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

! Utility functions in support of PIO io interface
module scream_scorpio_utils

  use iso_c_binding

  implicit none
  private
  save

  public :: scream_scorpio_createfile

  contains

  subroutine scream_scorpio_createfile()  bind(c)
    use pio,  only: pio_double 
!                    PIO_NOERR, PIO_CLOBBER

!    type(file_desc_t)  :: file
!    integer            :: ierr

    print *, 'I need some information...'
    print *, 'Just smile if you can hear me', pio_double

    ! Create new file
!    ierr = pio_createfile(pio_subsystem, file, pio_iotype,
!    "scream_pio_test", PIO_CLOBBER)

  end subroutine scream_scorpio_createfile

end module scream_scorpio_utils
