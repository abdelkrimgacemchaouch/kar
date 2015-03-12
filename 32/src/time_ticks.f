      subroutine time_ticks(count, count_rate, count_max)

      integer*4 count, count_rate, count_max

c Wall clock time in ticks.
c count is the number of ticks
c count_rate is the number of ticks per second
c count_max is the maximum value of the counter
      CALL system_clock(count, count_rate, count_max)

      return
      end
