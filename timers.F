      SUBROUTINE init_timers
      use timers_m
      IMPLICIT none
c
      INTEGER k
c
      do k=1,max_timers
        it_cnt(k) = 0
        tm_min(k) = 1.d12
        tm_max(k) = -1.d0
        tm_sum(k) = 0.d0
      enddo
c
      return
      END
c
      SUBROUTINE begin_timer(k)
      use timers_m
      IMPLICIT none
      INTEGER k
c
cc  rtc must be an 8-byte real
c      REAL*8 rtc
c      EXTERNAL rtc
c
      INTEGER*8 count
c
c      time1(k) = rtc()
c
      call system_clock(count)
      time1(k) = dble(count)
c
      return
      END
c
      SUBROUTINE end_timer(k)
      use timers_m
      IMPLICIT none

      INTEGER k
c
      REAL*8 tdiff
cc  rtc must be an 8-byte real
c      REAL*8 rtc
c      EXTERNAL rtc
c
      INTEGER*8 count,count_rate
c
      it_cnt(k) = it_cnt(k) + 1
c      tdiff = rtc() - time1(k)
      call system_clock(count,count_rate)
      tdiff = (dble(count) - time1(k))/dble(count_rate)
      if (tdiff.lt.tm_min(k)) tm_min(k) = tdiff
      if (tdiff.gt.tm_max(k)) tm_max(k) = tdiff
      tm_sum(k) = tm_sum(k) + tdiff
c
      return
      END
c
      SUBROUTINE print_timers
      use timers_m
      IMPLICIT none
c
      INTEGER k,n
c
      do k=1,max_timers
        n = it_cnt(k)
        if (n.gt.0) print 1000,
     &    k, n, tm_sum(k)/dble(n),tm_min(k),tm_max(k)
      enddo
c
      return
 1000 FORMAT('Timer',i4,':',i8,' calls, mean =',f10.5,
     &  ', min =',f10.5,
     &  ', max =',f10.5)
      END
