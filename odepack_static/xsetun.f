      subroutine xsetun (lun)
c
c this routine resets the logical unit number for messages.
c
      integer lun, mesflg, lunit
      common /eh0001/ mesflg, lunit
c
      !$OMP THREADPRIVATE(/eh0001/)
      if (lun .gt. 0) lunit = lun
      return
c----------------------- end of subroutine xsetun ----------------------
      end
