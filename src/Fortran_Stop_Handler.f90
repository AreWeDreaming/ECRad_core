! Dummy subroutine to prevent link error when
! projects for following console application are built:
! travisECRH, travisECE, travisOXB...
! The real routine needed for ECE web-service
! is specified inside ECE web-service project 
! using setjmp/longjmp exeption technique.
  SUBROUTINE FORTRAN_STOP_HANDLER(message)
   CHARACTER(Len=*), INTENT(in) :: message
   RETURN
  END
