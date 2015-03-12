      subroutine zonevols
c***********************************************************************
c     calculates absolute cell volumes (cc) needed for plankian source.                           *
c***********************************************************************
      include 'params.inc'
      include 'geomz.inc'
      
      do 20 j = 1,nzones
         volcl(j) = 0.0d0
 20   continue
      
      do 21 ks = 1,4
         do 22 j = 1,nzones
            if (itype(j,ks) .ne. 2 .and. itype(j,ks) .ne. 4) then
               tv1 = (zz(j,ks+1) - 
     .              zz(j,ks))*(rr(j,ks+1)*(rr(j,ks+1)
     .              + rr(j,ks)) + rr(j,ks)*rr(j,ks))
               volcl(j) = volcl(j) + 1.0471976*tv1
            endif
 22      continue
 21   continue
      
      return
      end
