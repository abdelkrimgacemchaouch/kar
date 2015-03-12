      subroutine copyseed( myIter, mySeed, firstRanf )

      include 'params.inc'
      include 'randseed.inc'

      integer myIter, myseed(*)
      REAL*4 firstRanf
      REAL*4 ranf
      external ranf

      do 1 i = 1, IrandNumSize
         mySeed(i) = seedarray(i,myIter)
    1 continue


      firstRanf = ranf(mySeed)

      return
      end
