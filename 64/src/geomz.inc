c This is 'geomz.inc'
c
c      scatter, shared
      dimension   
     .    incr (5),
     .    rr   (nzmax,5),
     .    zz   (nzmax,5),
     .    itype(nzmax,5),
     .    bom(nzmax,5),
     .    sqm(nzmax,5),
     .    ximp (nzmax),
     .    mid  (nzmax),
     .    volcl(nzmax),
     .    ximpl(nrzmax)

c      common, shared /geomz/ naxl,axl,nradl,radl,nzones
      common /geomz/ axl,naxl,nradl,radl,nzones
      

c      common, private/pgeomz/
      common /pgeomz/
     .    rr   ,
     .    zz   ,
     .    bom  ,
     .    sqm  ,
     .    ximp ,
     .    mid  ,
     .    volcl,
     .    ximpl,
     .    itype,
     .    incr

 
