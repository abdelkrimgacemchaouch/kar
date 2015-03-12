      program sphotmk
c***********************************************************************
c     sphotmk - scalar photon transport micro-kernel
c
c     this program is a scalar version of the vectorized monte carlo
c     code vphot which is used to solve photon transport problems in
c     spherical geometry.  this version INCLUDEs a plankian source,
c     thomson scattering, russian roulette/splitting and time census.
c     there are two options for russian roulette, one is based on photon
c     weight (culling) and the other is based on cell importances.
c***********************************************************************
c     * * * *  sphotmk - I/O units * * * *
c     
c     4  = 'input.dat'- Input File
c     6  = sdtout
c***********************************************************************
 
      INCLUDE 'params.inc'
      INCLUDE 'globals.inc'
      INCLUDE 'geomz.inc'
      INCLUDE 'randseed.inc'
      INCLUDE 'times.inc'

      INTEGER mySeed(IRandNumSize)
      INTEGER*4 nescgp(negrps)
      REAL*4 enesc(negrps)
      REAL*4 wmin, wmax
      REAL*4 wlost, wesc, wrr
      REAL*4 wabs, wcen, epgain, etot
      REAL*4 ranf
      EXTERNAL ranf

      REAL*4 t_enescSum, avg_enescSum, tmp_wesc

      INTEGER*4 nphtot, nploss
      INTEGER*4 nlost, nesc, nrr, nabs, ncen
      INTEGER*4 nscat, nsplt

      REAL*4  depArray(nzmax), depBuff(nzmax)
      REAL*4  r1, fRanBuf

      INTEGER*8 tick_count_start, tick_count_end
      INTEGER*8 tick_count_rate, tick_count_max

      CALL rdinput(nRuns)

c.....If the number of runs > maxruns then quit
      IF( nRuns .gt. maxruns ) THEN
         PRINT *, "nRuns > maxruns.  Need to modify maxruns in
     & params.inc file and recompile."
         PRINT *, "exiting..."
         goto 1001
      END IF

      CALL seedranf(nRuns+1)
      CALL genmesh
      CALL genxsec

      DO i = 1, nzones
         depArray(i) = 0.D0
         depBuff(i)  = 0.D0
      END do

      t_enescSum = 0.0
      ithRun = 1

      CALL time_ticks(tick_count_start, tick_count_rate, tick_count_max)

c.....main loop
      DO 1000 ict = 1, nRuns

         CALL execute(ithRun, mySeed, nescgp, enesc, 
     &        wmin, wmax, wlost, wesc, wrr, wabs, wcen, epgain, etot,
     &        nphtot, nploss, nlost, nesc, nrr, nabs, ncen,
     &        nscat, nsplt, fRanBuf, depArray )

         tmp_wesc = 0.0
         do 498 i=1,12
             tmp_wesc = tmp_wesc + enesc (i)
 498     continue

         t_enescSum = t_enescSum + ( tmp_wesc/dble(nphtot) )
         r1 = ranf(mySeed)
         ithRun = ithRun + 1

 1000 CONTINUE

      CALL time_ticks(tick_count_end, tick_count_rate, tick_count_max)

      avg_enescSum = t_enescSum/nRuns

      tick_count = tick_count_end - tick_count_start
      if(tick_count < 0) tick_count = tick_count + tick_count_max
      progWallClockTime = tick_count / tick_count_rate


c.....Write results:
c
c.....number of iterations       -> set in rdinput
c.....average escape probability -> final result to check for correctness
c.....execution time             -> time in sec (from time_ticks.f)
c
c.....Here the results are written to stdio but
c.....this can be modified according to the environment.
c.....For example for hardware simulator these can be stored
c.....in registers instead of being printed out.
       write(6,600) NRuns, avg_enescSum, progWallClockTime
 600   format('Sequoia Benchmark Version 1.0',    
     &        //, 'NRuns           = ', i20, 
     &        //, 'avg. esc. prob. = ', f20.6, 
     &        //, 'seconds         =', f20.2)


 1001 CONTINUE

      STOP
      END
