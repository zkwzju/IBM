      module timers_m
      INTEGER max_timers
      PARAMETER (max_timers=200)
      INTEGER it_cnt(max_timers)
      REAL*8 tm_sum(max_timers),tm_min(max_timers),tm_max(max_timers)
c time1 must be 8-byte real array
      REAL*8 time1(max_timers)
      COMMON /timers/ it_cnt,tm_sum,tm_min,tm_max,time1
      end module
