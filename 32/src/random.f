c==============================================================================
c This file contains all the subroutines and functions related to
c random number generation
c==============================================================================


c==============================================================================
c------------------------------------------------------------------------------
    	 real*4 function ranf( mySeed )
c------------------------------------------------------------------------------
      include 'params.inc'
      include 'randseed.inc'
      INTEGER mySeed(IRandNumSize)
      real*4  randnum
      CALL pranf(mySeed,randnum)
      ranf = randnum
      return
      end


c==============================================================================
c------------------------------------------------------------------------------
      subroutine seedranf( numStreams )
c------------------------------------------------------------------------------
      include 'randseed.inc'
      integer*4 numStreams
      integer*4 tempNumStreams  !Because rans may change numStreams.

      NStreams = numStreams
      tempNumStreams = numStreams

      CALL rans(tempNumStreams,0,seedarray)

      return 
      end


c==============================================================================
c------------------------------------------------------------------------------
      subroutine pranf( Seed, RandNum )
c------------------------------------------------------------------------------
c     
c     This is the one of the two user-callable routines of a linear congruential
c     random number generator, modeled after the RANF generator for the Cray.
c     This routine generates the next random number in the sequence and a new seed
c     for the remainder of the sequence.  The seed and the random number are the
c     same, but are returned in different form: the random number is a fortran
c     'real', but the seed is an array of four words, each containing an integer
c     that is used internally to the generator as one digit of a four-digit,
c     modulo-4096 integer.
c     
c     It returns the new random number and a new seed.
c     
c------------------------------------------------------------------------------
      include 'pranf.inc'
      integer Seed( IRandNumSize )
      real*4 RandNum

c**** Data common to the PRanf package.
      integer Multiplier( IRandNumSize ), DefaultSeed( IRandNumSize )
      real*4 Divisor( IRandNumSize )

      common / PRanfCom / Multiplier, DefaultSeed, Divisor 


c**** End of PRanf common data

      RandNum = float( Seed( 4 ) ) / Divisor( 4 ) +
     1     float( Seed( 3 ) ) / Divisor( 3 ) +
     2     float( Seed( 2 ) ) / Divisor( 2 ) +
     3     float( Seed( 1 ) ) / Divisor( 1 ) 

      CALL ranfmodmult( Multiplier, Seed, Seed )

      return
      end

      BLOCK DATA PINIT

      include 'pranf.inc'
      integer Multiplier( IRandNumSize ), DefaultSeed( IRandNumSize )
      real*4 Divisor( IRandNumSize )

      common / PRanfCom / Multiplier, DefaultSeed, Divisor 

      data Multiplier / 373, 3707, 1442, 647 /
      data DefaultSeed / 3281, 4041, 595, 2376 / 
      data Divisor / 281474976710656.0, 68719476736.0, 16777216.0, 
     &                  4096.0 /

      END


c==============================================================================
c------------------------------------------------------------------------------
      subroutine Rans( NIn, Seed1, SeedArray )
c------------------------------------------------------------------------------
c     
c     This routine divides the sequence of random numbers into N subsequences,
c     each with its own seed.  The seeds for the independent subsequences are
c     returned in the seed array.  if Seed1 is zero, all zeroes will be returned.
c     To prevent this, Seed1 is set to [3281, 4041, 595, 2376], which is 
c     statistically the best starting seed.  The wheel is then divided into the
c     N pieces (where N is odd and >= NIn) by dividing its period (2**46) by N.
c     
c     Then, seed(i) = seed(i-1) * (a**k mod 2**48), and 1<=k<=N.
c     
c     Here, 'a' is the multiplier used by the linear congruential generator whose
c     wheel we are dividing up.
c     
c     The number of streams must be odd; if NIn is even N will be NIn+1, and
c     n extra stream of random numbers will be available that will not get used.
c     
c     It returns an array of seeds, each an array of 4 integers that are used as
c     the digits of a four-digit modulo-4096 integer.
c     
c------------------------------------------------------------------------------
      include 'pranf.inc'

      integer NIn, Seed1, SeedArray, N, atothek, K, 
     1     KBinary, InSeed, OutSeed

c**** Data common to the PRanf package.
      integer Multiplier, DefaultSeed

      real*4 Divisor

      dimension Multiplier( IRandNumSize ), DefaultSeed( IRandNumSize ),
     1     Divisor( IRandNumSize )

      common / PRanfCom / Multiplier, DefaultSeed, Divisor 
c**** End of PRanf common data

      dimension SeedArray( IRandNumSize, MaxStreams ),
     1     atothek( IRandNumSize ),
     2     K( IRandNumSize ),
     3     KBinary( IBinarySize ),
     4     InSeed( IRandNumSize ),
     5     OutSeed( IRandNumSize )

      INTEGER NisOdd

c     scatter, shared SeedArray
c------------------------------------------------------------------------------

      IZero = 0

c**** Make sure we are generating an odd number of random number sequences.
      if( NisOdd( NIn ) .eq. 1 ) then
         N = NIN
      else
         N = NIn + 1
      endif

c**** Set up the initial seed to either a legal input value or its default
c**** values. The input integer, if nonzero, is used for the first of the
c**** four modulo-4096 digits of the actual initial seed.

c**** The First element of the returned SeedArray will be used here, since
c**** at least one seed will be returned, and the first seed of the set 
c**** returned will be the seed for the entire wheel.
      if( Seed1 .eq. IZero ) then
         SeedArray( 1, 1 ) = DefaultSeed( 1 )
         SeedArray( 2, 1 ) = DefaultSeed( 2 )
         SeedArray( 3, 1 ) = DefaultSeed( 3 )
         SeedArray( 4, 1 ) = DefaultSeed( 4 )
      else
         SeedArray( 1, 1 ) = abs( Seed1 )
         SeedArray( 2, 1 ) = IZero
         SeedArray( 3, 1 ) = IZero
         SeedArray( 4, 1 ) = IZero
      endif

c**** 'a' is the multiplier for the Ranf linear congruential generator.
      if( N .eq. 1 ) then
c*******If only one stream of random numbers is needed, do not bother to 
c*******raise 'a' to the first power.
         atothek( 1 ) = Multiplier( 1 )
         atothek( 2 ) = Multiplier( 2 )
         atothek( 3 ) = Multiplier( 3 )
         atothek( 4 ) = Multiplier( 4 )
      else
c*******more than one stream is needed; generate the Kth seed by multiplying
c*******the K-1st seed by the multiplier raised to the Nth power.
         CALL ranfk( N, K )
         CALL ranfkbinary( K, KBinary )
         CALL ranfatok( Multiplier, KBinary, atothek )
         do 100 I = 2, N
            InSeed( 1 ) = SeedArray( 1, I-1 )
            InSeed( 2 ) = SeedArray( 2, I-1 )
            InSeed( 3 ) = SeedArray( 3, I-1 )
            InSeed( 4 ) = SeedArray( 4, I-1 )
            CALL ranfmodmult( InSeed, atothek, OutSeed )
            SeedArray( 1, I ) = OutSeed( 1 )
            SeedArray( 2, I ) = OutSeed( 2 )
            SeedArray( 3, I ) = OutSeed( 3 )
            SeedArray( 4, I ) = OutSeed( 4 )
 100     continue
      endif

 110  format( '   leaving Rans' )

      return
      end


c MJC mod to remove statement function
      INTEGER FUNCTION NISODD( M )
c****  check for argument numbers being odd.
      NIsOdd = mod( m, 2 )
      RETURN
      END


c==============================================================================
c------------------------------------------------------------------------------
      subroutine ranfatok( a, Kbinary, atothek )
c------------------------------------------------------------------------------
c     
c     This routine calculates a to the Kth power, mod 2**48. K is a binary number.
c     
c     It returns the calculated value as an array of four modulo-4096 digits.
c     
c------------------------------------------------------------------------------
      include 'pranf.inc'
      integer a, KBinary, atothek, asubi
      dimension a( IRandNumSize ),
     1     KBinary( IBinarySize ),
     2     atothek( IRandNumSize ),
     3     asubi( IRandNumSize )
c...MJC mod...this was missing
      INTEGER IZero
      PARAMETER( IZero = 0 )

c------------------------------------------------------------------------------
c**** The following amounts to the first iteration of a 46-loop.
      asubi( 1 ) = a( 1 )
      asubi( 2 ) = a( 2 )
      asubi( 3 ) = a( 3 )
      asubi( 4 ) = a( 4 )

      atothek( 1 ) = 1
      atothek( 2 ) = IZero
      atothek( 3 ) = IZero
      atothek( 4 ) = IZero

      do 100 I = 1, 45
         if( KBinary( I ) .ne. IZero ) then
            CALL ranfmodmult( atothek, asubi, atothek )
         endif
         CALL ranfmodmult( asubi, asubi, asubi )
 100  continue

      return
      end


c==============================================================================
c------------------------------------------------------------------------------
      subroutine ranfk( N, K )
c------------------------------------------------------------------------------
c     
c     This routine calculates 2**46/N, which should be the period of each of the
c     subsequences of random numbers that are being created. Both the numerator
c     and the result of this calculation are represented as an array of four
c     integers, each of which is one digit of a four-digit moduo-4096 number.  The
c     numerator is represented as (1024, 0, 0, 0 ), using base ten digits.
c     
c     It returns the result of the division.
c     
c------------------------------------------------------------------------------

      include 'pranf.inc'

      integer N, K, nn, r4, r3, r2, q4, q3, q2, q1

      dimension K( IRandNumSize )

c------------------------------------------------------------------------------
      nn = N + iranfeven( N )

      q4 = 1024 / nn
      r4 = 1024 - (nn * q4)
      q3 = (r4 * 4096) / nn
      r3 = (r4 * 4096) - (nn * q3)
      q2 = (r3 * 4096) / nn
      r2 = (r3 * 4096) - (nn * q2)
      q1 = (r2 * 4096) / nn

      K( 1 ) = q1
      K( 2 ) = q2
      K( 3 ) = q3
      K( 4 ) = q4

      return
      end


c==============================================================================
c------------------------------------------------------------------------------
      subroutine ranfkbinary( K, KBinary )
c------------------------------------------------------------------------------
c     
c     This routine calculates the binary expansion of the argument K, which is a
c     48-bit integer represented as an array of four 12-bit integers.
c     
c     It returns an array of 48 binary values.
c     
c------------------------------------------------------------------------------
      include 'pranf.inc'
      integer K, KBinary, X, Bits
      dimension K( IRandNumSize ),
     1     KBinary( IBinarySize ),
     2     Bits( Mod4096DigitSize )

c------------------------------------------------------------------------------
      do 300 I = 1, 4
         X = K( I ) / 2
         Bits( 1 ) = iranfodd( K( I ) )

         do 100 J = 2, Mod4096DigitSize 
            Bits( J ) = iranfodd( X )
            X = X / 2
 100     continue

         do 200 J = 1, Mod4096DigitSize
            KBinary( (I-1)*Mod4096DigitSize + J ) = Bits( J )
 200     continue

 300  continue

      return
      end


c==============================================================================
c------------------------------------------------------------------------------
      subroutine ranfmodmult( A, B, C )
c------------------------------------------------------------------------------
c     
c     Ths routine calculates the product of the first two arguments.  All three
c     arguments are represented as arrays of 12-bit integers, each making up
c     the digits of one radix-4096 number.  The multiplication is done 
c     piecemeal.
c     
c     It returns the product in the third argument.
c     
c------------------------------------------------------------------------------
      include 'pranf.inc'
      integer A, B, C, j1, j2, j3, j4, k1, k2, k3, k4
      integer a1,a2,a3,b1,b2,b3
      dimension A( IRandNumSize ),
     1     B( IRandNumSize ),
     2     C( IRandNumSize )

c------------------------------------------------------------------------------
c     j1 = A( 1 ) * B( 1 )
c     j2 = A( 1 ) * B( 2 ) + A( 2 ) * B( 1 )
c     j3 = A( 1 ) * B( 3 ) + A( 2 ) * B( 2 ) + A( 3 ) * B( 1 )
c     j4= A( 1 ) * B( 4 ) + A( 2 ) * B( 3 ) + A( 3 ) * B( 2 ) + A( 4 ) * B( 1 )
      a1 = A(1)
      a2 = A(2)
      a3 = A(3)
      b1 = B(1)
      b2 = B(2)
      b3 = B(3)
      j1 = a1 * b1
      j2 = a1 * b2 + a2 * b1
      j3 = a1 * b3 + a2 * b2 + a3 * b1
      j4 = a1 * B( 4 ) + a2 * b3 + a3 * b2 + A( 4 ) * b1

      k1 = j1
      k2 = j2 + k1 / 4096
      k3 = j3 + k2 / 4096
      k4 = j4 + k3 / 4096

      C( 1 ) = mod( k1, 4096 )
      C( 2 ) = mod( k2, 4096 )
      C( 3 ) = mod( k3, 4096 )
      C( 4 ) = mod( k4, 4096 )

      return
      end


c==============================================================================
c------------------------------------------------------------------------------
      function iranfeven( N )
c------------------------------------------------------------------------------
c     
c     This function checks the parity of the argument integer.
c     
c     It returns one if the argument is even and zero otherwise.
c     
c------------------------------------------------------------------------------
      include 'pranf.inc'
      integer N

c------------------------------------------------------------------------------
      if( mod( N, 2 ) .eq. 0 ) then
         iranfeven = 1
      else
         iranfeven = 0
      endif

      return
      end


c==============================================================================
c------------------------------------------------------------------------------
      function iranfodd( N )
c------------------------------------------------------------------------------
c     
c     This function checks the parity of the argument integer.
c     
c     It returns one if the argument is odd and zero otherwise.
c     
c------------------------------------------------------------------------------
      include 'pranf.inc'
      integer N

c------------------------------------------------------------------------------
      if( mod( N, 2 ) .eq. 1 ) then
         iranfeven = 0
      else
         iranfeven = 1
      endif

      iranfodd = 1 - iranfeven
      return
      end
