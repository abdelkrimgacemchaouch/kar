c------------------------------------------------------------------------------
c
c All the input parameters are set in this routine. If a parameter
c is changed here this routine must be recompiled and re-linked.
c The input data is not read from an input file in order to eliminate
c file I/O and therefore make it faster and easier to use in hardware 
c simulation environments
c
c------------------------------------------------------------------------------
      subroutine rdinput( NRuns )

      include 'params.inc'
      include 'geomz.inc'
      include 'globals.inc'

      INTEGER Nruns, matb1, mate1, matb2, mate2
      INTEGER nout
      parameter (nout=10)


c------------------------------------------------------------------------------
c     NRuns and bwgt affect the program execution time directly.
c     Edit these as needed
c------------------------------------------------------------------------------

c---- The number of iterations in the outermost loop.
c---- Execution time increases linearly with NRuns
      NRuns = 4

c---- The bundle weight (kev).
c---- The time per iteration decreases with increasing bwgt
      bwgt = 3.12d+14
c------------------------------------------------------------------------------


c------------------------------------------------------------------------------
c     The following parameters should not be modified
c------------------------------------------------------------------------------

c---- library usage
      ilib = 1
      illib = 0

c---- # of particles (different if plankian used)
      npart = 4000

c---- energy bins (0=12, 1-12=1, 13=ross.mean)
      igroup = 0

c---- opacity (must be set to 2)
      ixopec = 2

c---- source (1=uniform in sphere, 2=plankian)
      isorst = 0

c---- r-roulette/split (0=none, 1=impt, 2=size)
      irr = 0

c---- thomson scattering (0=not used, 1=used)
      ithom = 0  

c---- number of axial meshes
      naxl = 49

c---- number of radial meshes
      nradl = 40

c---- number of material regions
      nreg = 2

c---- low weight cutoff
      wcut = 7.5d-01  

c---- time to census (sec)
      tcen = 1.0d-12  

c---- weight mult. for russian roulette
      xmult = 1.05d+00
      
c---- portion of sphere analyzed (degrees)'
      axl = 1.8d+02

c---- sphere radius (cm)'
      radl = 1.0d-03

c---- input opacity (1/cm)'
      opec = 5.0d-01


c---- For the purposes of the ASCI Purple Benchmark, 
c---- it is mandatory that nreg = 2
      nreg = 2

c---- Set mid array
      matb1 = 1
      mate1 = 980
      matb2 = 981
      mate2 = 1960

      do 151 j = matb1, mate1
         mid(j) = 1
 151  continue
      do 152 j = matb2, mate2
         mid(j) = 2
 152  continue

c---- region 1 material
      mtl(1) = 3
c---- region 1 atomic ratio
      atrat(1) = 4.0d-01
c---- region 1 density(g/cc)
      dns(1) = 1.0d+03
c---- region 1 temperature(ev)
      tmp(1) = 2.0d+04

c---- region 2 material
      mtl(2) = 2
c---- region 2 atomic ratio
      atrat(2) = 5.0d-01
c---- region 2 density(g/cc)
      dns(2) = 1.0d+02
c---- region 2 temperature(ev)
      tmp(2) = 1.0d+03
      
      return
      end
