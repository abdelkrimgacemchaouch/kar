      subroutine genxsec
c----------------------------------------------------------------------
c     Setup the opacity data by group (13) and region (nreg).
c     Use opacities from data statement (sigdat).
c----------------------------------------------------------------------
      include 'params.inc'
      include 'geomz.inc'
      include 'globals.inc'

      REAL *4 sigdat(2,negp1)

      data sigdat/     .1395449d+06,.8093366d+06,.2022023d+05,
     |     .1344934d+06,.4188673d+04,.4467951d+05,.1033437d+04,
     |     .3250243d+05,.4297604d+03,.1237669d+05,.2235549d+03,
     |     .5984154d+04,.1270899d+03,.3111531d+04,.7992034d+02,
     |     .1722049d+04,.4425178d+02,.9100233d+03,.1456625d+02,
     |     .3776509d+03,.7138723d+00,.7534896d+02,.4231184d-01,
     |     .1032693d+01,.7666376d+01,.6094672d+04/
      
      xthom = 0.0d0
c-----ithom non-zero implies that thomson scattering is used
      if (ithom .ne. 0) xthom = 1.0d0
      
      do 75 i = 1,nreg         
c-----store the thomson scattering cross section.
c     - thomson cross section depends only on the material,
c     - and not on the energy group         
         sigth(i) = .4006*atrat(i)*dns(i)*xthom

         do 69 ig = 1,13            
c-----store total cross section without thomson.
c     - first index in cross section is the material id,
c     - the second is the energy group
            sigtot(i,ig) = sigdat(i,ig)            
            sig  (i,ig) = sigtot(i,ig) + sigth(i)
            scrat(i,ig) = sigth(i) / sig(i,ig)
 69      continue
 75   continue

      return
      end
