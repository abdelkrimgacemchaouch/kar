c This is 'globals.inc'
c
c      scatter, shared 
       dimension
     .    mtl  (nmrmax),
     .    dns  (nmrmax),
     .    tmp  (nmrmax),
     .    atrat(nmrmax),
     .    sig   (nmrmax,negp1),
     .    scrat (nmrmax,negp1),
     .    sigtot(nmrmax,negp1),
     .    sigth (nmrmax),
     .    opcv  (       negp1)

c      common, shared /material/ nreg
      common  /material/ nreg

c      common, private /pmaterial/
      common /pmaterial/
     .    mtl  ,
     .    dns  ,
     .    tmp  ,
     .    atrat
c      common, shared /xsec/ opec
      common /xsec/ opec

c      common, private /pxsec/
      common /pxsec/
     .    sig   ,
     .    scrat ,
     .    sigtot,
     .    sigth ,
     .    opcv  
 
c      common, shared /controls/ igroup,ixopec,isorst,irr,ithom,
      common /controls1/ igroup,ixopec,isorst,irr,ithom,ilib,illib
      common /controls2/ tcen,xmult,title(20)
 
c      common, shared /newcom/ bwgt    !Formerly 1st item in /samples/.
      common /newcom/ bwgt    !Formerly 1st item in /samples/.


c      common /samples/ nescgp(negrps), enesc(negrps), wcut, wmin,  wmax,
c     .   wlost, wesc, wrr, wabs, wcen, epgain, etot,
c     .   npart, nphtot, nploss, 
c     .   nlost, nesc, nrr, nabs, ncen,
c     .   nscat, nsplt, ntrac


