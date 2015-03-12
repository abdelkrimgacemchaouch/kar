      subroutine genmesh
c----------------------------------------------------------------------
c     A simple polar mesh is generated.
c----------------------------------------------------------------------
      include 'params.inc'
      include 'geomz.inc'
      include 'globals.inc'
      
      real*4 m
c--------------------initialize variables
      pi      = 3.141592653589793238d0
      incr(1) = 1
      incr(2) = -naxl
      incr(3) = -1
      incr(4) = naxl
      
      if (irr .ne. 0) then
c--------------------assign cell importances (or weight limits) for rr
       
cc.....Seems that impl decreases geometrically with radius....MJC
         if (irr .eq. 1) then

            ximpl(1) = 1.0d0
            do 221 i = 2,nradl
               ximpl(i) = xmult * ximpl(i-1)
 221        continue
            
c...fill a polar array with these values, equal azumuthally MJC
            do 222 i = 1,nradl
               do 223 k = 1,naxl
                  iz       = k + (i-1)*naxl
                  ximp(iz) = ximpl(i)
 223           continue
 222        continue
            
         else
            wmin = wcut
            wmax= 1.0d0 / wcut
         endif
         
      endif

c..The number above implies that the center node is #1 and that the next
c..layer of nodes around is is numbered increasingly either clockwise or
c..counter clockwise. MJC

      
c--------------------assign cell numbers and cell coordinates
      
      nzones = nradl * naxl
      hr     = radl  / nradl
      ha     = (axl*pi)/(naxl*180.0d0)
      naxpls = naxl + 1
      do 1 iang = 1,naxpls
         ang = (iang-1)*ha
         cs  = cos(ang)
         sn  = sin(ang)
         do 2 ird=1,nradl
            r = (ird*hr)*sn
            z = (ird*hr)*cs
            
            if (iang .le. 1) goto 50
            jc       = ird*naxl - iang + 2
            zz(jc,4) = z
            rr(jc,4) = r
            
            if (ird .ge. nradl) goto 50
            jc       = jc+naxl
            zz(jc,3) = z
            rr(jc,3) = r
            
 50         continue
            
            if (iang .ge. naxpls) goto 51
            jc       = ird*naxl - iang + 1
            zz(jc,1) = z
            rr(jc,1) = r
            if (ird .ge. nradl) goto 51
            jc       = jc+naxl
            zz(jc,2) = z
            rr(jc,2) = r
            
 51         continue
            
 2       continue
 1    continue
      
      do 3 ii = 1,nzones
         do 8 ij = 1,4
            if (abs(zz(ii,ij)) .lt. 1.0e-9) zz(ii,ij) = 0.0d0
            if (abs(rr(ii,ij)) .lt. 1.0e-9) rr(ii,ij) = 0.0d0
 8       continue
         zz(ii,5) = zz(ii,1)
         rr(ii,5) = rr(ii,1)
 3    continue


c------------------------------compute SLOPE and INTERCEPT of cell sides,
c     and categorize special cases:
c     1 = normal cell,
c     2 = cell on z-axis,
c     3 = horizontal side,
c     4 = vertical side
      
      do 60 i = 1,nzones
         do 61 ks = 1,4
            rdif = rr(i,ks+1) - rr(i,ks)
            zdif = zz(i,ks+1) - zz(i,ks)
            
            if ((rr(i,ks).eq.0.0d0).and.(rr(i,ks+1).eq.0.0d0)) then
               itype(i,ks) = 2
            else if (abs(rdif) .lt. 1.0d-9) then
               itype(i,ks) = 3
            else if (abs(zdif) .lt. 1.0d-9) then
               itype(i,ks) = 4
            else
               itype(i,ks) = 1
            end if
            
            if (itype(i,ks) .eq. 1) then
               
               if ((ks .eq. 1) .or. (ks .eq. 3)) then
                  m = rdif/zdif
                  b = 0.0d0
                  sqm(i,ks) = m*m
                  bom(i,ks) = 0.0d0
                  
               else
                  b = (rr(i,ks)*zz(i,ks+1)
     &                 -rr(i,ks+1)*zz(i,ks))/zdif
                  m = rdif/zdif
                  sqm(i,ks) = m*m
                  bom(i,ks) = b/m
               endif
               
            endif
 61      continue
 60   continue
      
      CALL zonevols
      
      return
      end
