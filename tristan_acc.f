C===================== beginning of tristan.f ==========================
C  -->  In this version of TRISTAN particles which hit the boundaries
C  -->  are arrested there and redistributed more uniformily by having
C  -->  the boundaries slightly conducting: this allows electrons to
C  -->  recombine with ions and provides a realistic way of eliminating
C  -->  escaping particles from the code. Fresh particle fluxes can
C  -->  then be introduced independently across the boundaries.
C
C
C
C  TRIdimensional STANford code, TRISTAN, fully electromagnetic,
C  with full relativistic particle dynamics. Written during spring
C  1990 by OSCAR BUNEMAN, with help from TORSTEN NEUBERT and
C  KEN NISHIKAWA, and with support from a NASA contract.
      real me,mi
      logical in
      common /partls/ x(110000),y(110000),z(110000),
     &u(110000),v(110000),w(110000)
      common /fields/ ex(65,35,35),ey(65,35,35),ez(65,35,35),
     &bx(65,35,35),by(65,35,35),bz(65,35,35),c,ix,iy,iz
C  The sizes of the particle and field arrays must be chosen in
C  accordance with the requirements of the problem and the limitation
C  of the computer memory. Any changes from the choices presented 
C  here must also be applied in the COMMON statements of the "mover"
C  and "depsit" subroutines. There the field arrays are treated as
C  single-indexed of length "lot"=65*35*35=79625.
      common /spread/q,sm(27),ms(27)
      dimension table(64)
      write(*,'(6h last?)')
      read(*,'(i4)')last
      rewind(9)
      if(last.gt.64)go to 5
C  If "last" (the last timestep) is less or equal 64, the following
C  initialisations are performed prior to time-stepping.
      open(unit=7,file='flds.d',status='new',form='unformatted')
      open(unit=8,file='prtl.d',status='new',form='unformatted')
C  The ions occupy the first half, the electrons the second half of
C  the particle arrays whose total length is:
      maxptl=110000
      maxhlf=maxptl/2
C  The dimensions of the field arrays are:
      mx=65
      my=35
      mz=35
C  For CRAY-s, the first two field dimensions should be ODD.
C  Fields can be treated as triple-indexed or single-indexed. In the
C  latter case, Fortran's "strides" for the three dimensions are:
      ix=1
      iy=mx
      iz=iy*my
      lot=iz*mz
C  Initialise some physical constants:
      qe=-1.
      qi=1.
      me=1.
      mi=16.
      qme=qe/me
      qmi=qi/mi
      c=.5
C Our finite difference equations imply delta_t = delta_x =
C delta_y = delta_z = 1. So c must satisfy the Courant condition.
C Our bx,by and bz arrays are really records of c*Bx, c*By,
C c*Bz: this makes for e <-----> b symmetry in Maxwell's equations.
C Further, field units are such that epsilon_0  is 1.0 and hence mu_0
C is 1/c**2. This means that for 0.8 electrons per cell (see example
C of particle initialisation below) omega_p-squared is 0.8*qe**2/me.
C ---------------
C  Initialise the fields, preferably both E and B with zero divergence
C  so that Poisson's equation (only an initial condition in our charge
C  conserving code) is easily satisfied by initialising ions and 
C  electrons in the same locations. A common initialisation of the
C  fields is to uniform vectors, such as a uniform bz or a uniform
C  bz with a uniform ey (as seen in a frame of reference moving 
C  x-ward across bz). Remember that our bx,by,bz are really cBx,cBy,
C  cBz and omega-c for the electrons and ions is qme|b|/c, qmi|b|/c.
C -----------
C  In this example we initialise all fields to zero and let Earth's
C  magnetic field get established transiently by a current loop placed
C  at the center of the Earth (see loop 66):
      do 1 k=1,mz
      do 1 j=1,my
      do 1 i=1,mx
      ex(i,j,k)=0.
      ey(i,j,k)=0.
      ez(i,j,k)=0.
      bx(i,j,k)=0.
      by(i,j,k)=0.
   1  bz(i,j,k)=0.
C  The b-components can, if desired, be initialised to an extra-
C  terrestrial magnetic field. This will be traveling with the
C  solar wind and give rise to an apparent electric field which
C  must be initialised in the e-components.
C  The ring current "o" determines the strength of the Earth's dipole
C  field. It is made to evolve from zero to its full value within the
C  first 64 timesteps, smoothly and gently, to avoid large transients.
C  This is achieved by building up "o" from its first, second and third
C  time differences, o1,o2,o3, following a cubic. Initially:
      o3=-.25
      o2=-31.*o3
      o1=o2
      o=o1
C  The ratio -o2/o3 should be an exact integer n (here 31). The number
C  of time steps taken for the build-up is 2*n. The final value of "o"
C  will be -o3*n*(n+1)*(2*n+1)/3 .
C  The location of the Earth's dipole (along z, perp to Sun-Earth
C  axis) is at x=ie+.5,y=je+.5,z=ke where:
      ie=20
      je=17
      ke=18
C -----------
C Data for smoothing: the currents fed into Maxwell's equations
C are smoothed by convolving with the sequence  .25, .5, .25 in each
C dimension. Generate the 27 weights ("sm") and index displacements
C ("ms") which use the Fortran multiple index convention:
      n=1
      do 4   k=-1,1
      do 4   j=-1,1
      do 4   i=-1,1
      sm(n)=.015625*(2-i*i)*(2-j*j)*(2-k*k)
      ms(n)=ix*i+iy*j+iz*k
  4   n=n+1
C -----------
C  Many computers call for initialisation of their random number
C  generator with a "lucky number" seed:
      luck=12345
C  In the following, ranf() is a random real betweem 0. and 1.0. Use an
C  editor to replace all occurrences of "ranf()" by the appropriate random
C  number call of your system; 'ranf()' is the call for CRAYs. On DEC
C  stations we have used 'ran(luck)' instead.
C  The following table is specific to the solar wind problem: it is
C  used for picking particle velocities in the transverse solar
C  fluxes at the lateral boundaries.
      do 3   n=1,64
  3   table(n)=.5*sqrt(log(16384.)-2.*log(float(2*n-1)))
C --------------
C  Particle initialisation: The domain in which particles are allowed
C  is given by 3.0 =< x < mx-2, 3.0 =< y < my-2, 3.0 =< z < mz-2  .
C  Field indices are ifix(x), ifix(y) and ifix(z); this leaves three
C  empty guard cells at the top, two at the bottom in each dimension.
C  (Guard cells are needed because the radiation absorbing boundary
C  condition of subroutine SURFACE assumes vacuum in the outermost
C  field cells and because the current deposit of subroutine DEPSIT
C  spills over into neighboring cells, partly due to the smoothing
C  algorithm).
C     Here we initialise the solar wind to uniform density of 0.8
C  ion-electron pairs per cell (note spacing of 1.25 units along x)
C  and to velocities given by:
      vdrft=.25
C  (this will be in the x-direction)
      vth2e=.1
      vth2i=.025
C  (these are twice the electron- and ion- thermal velocities in the
C  sense of rms deviations of each single component)
      n=0
      do  2 x0=3.,mx-3.25,1.25
      do  2 y0=3.,my-3.,1.
      do  2 z0=3.,mz-3.,1.
      n=n+1
C  Add some random spread to these regular spacings:
      x(n)=x0+1.25*rand()
      y(n)=y0+rand()
      z(n)=z0+rand()
C  Place electrons in the same locations as ions for zero charge
C  density (consistent with zero or uniform inital electric fields):
      x(maxhlf+n)=x(n)
      y(maxhlf+n)=y(n)
      z(maxhlf+n)=z(n)
C  Adding three random numbers, each uniformily random between -.5
C  and +.5, leads to a distribution closely resembling maxwellian.
C  The factor in front then is twice the rms value of a single
C  "thermal" velocity component.
      u(n)=vdrft+vth2i*(vth2i+rand()+rand()+rand()-1.5)
      v(n)=vth2i*(rand()+rand()+rand()-1.5)
      w(n)=vth2i*(rand()+rand()+rand()-1.5)
      u(maxhlf+n)=vdrft+vth2e*(vth2e+rand()+rand()+rand()-1.5)
      v(maxhlf+n)=vth2e*(rand()+rand()+rand()-1.5)
   2  w(maxhlf+n)=vth2e*(rand()+rand()+rand()-1.5)
C  Care should be taken when chosing drift and thermal velocities
C  that none of the resultant velocities exceed c.
C  Record the total number of ions ("ions") and electrons ("lecs"):
      ions=n
      lecs=n
      lap=0
C  After these initialisations, enter the timestep loop at the point
C  where Earth's dipole field is established:
      go to 60
C  When continuing a run, begin by reading in the most recent field
C  and particle data on record:
   5  open(unit=7,file='flds.d',status='old',form='unformatted')
      open(unit=8,file='prtl.d',status='old',form='unformatted')
      read(7)ex,ey,ez,bx,by,bz,sm,ms,c,mx,my,mz,ix,iy,iz,lot,
     &o,o1,o2,o3,ie,je,ke,qe,qi,qme,qmi,vdrft,vth2e,vth2i,table,luck,lap
      read(8)ions,lecs,maxptl,maxhlf,x,y,z,u,v,w
      rewind(7)
      rewind(8)
C  The time-stepping loop:
   6  lap=lap+1
C Before moving particles, the magnetic field is Maxwell-advanced
C by half a timestep (after a preliminary partial update of the edges):
      call preledge(by,bz,bx,ey,ez,ex,iy,iz,ix,my,mz,mx,1,c)
      call preledge(bz,bx,by,ez,ex,ey,iz,ix,iy,mz,mx,my,1,c)
      call preledge(bx,by,bz,ex,ey,ez,ix,iy,iz,mx,my,mz,1,c)
      do 7 i=1,mx-1
      do 7 j=1,my-1
      do 7 k=1,mz-1
      bx(i,j,k)=bx(i,j,k) + (.5*c) *
     &          (ey(i,j,k+1)-ey(i,j,k)-ez(i,j+1,k)+ez(i,j,k))
      by(i,j,k)=by(i,j,k) + (.5*c) *
     &          (ez(i+1,j,k)-ez(i,j,k)-ex(i,j,k+1)+ex(i,j,k))
  7   bz(i,j,k)=bz(i,j,k) + (.5*c) *
     &          (ex(i,j+1,k)-ex(i,j,k)-ey(i+1,j,k)+ey(i,j,k))
C Now move ions and electrons:
      call mover(1,ions,qmi)
      call mover(1+maxhlf,lecs+maxhlf,qme)
C The Maxwell-advance of the fields begins with another half-step
C advance of the magnetic field since for Maxwell's equations the
C B - information and the E - information must be staggered in time.
C In space, their information is also staggered. Here we show the
C locations where field components are recorded:
C
C
C
C                z
C                ^
C                |
C                |
C                |
C                *---------------Ey---------------*
C               /|                                |
C              / |                                |
C             /  |                                |
C            Ex  |                                |
C           /    |                                |
C          /     |                                |
C         /      |                                |
C        *       |                                |
C        |       |                                |
C        |       |                                |
C        |       |                                |
C        |       Ez              Bx               Ez
C        |       |                                |
C        |       |                                |
C        |       |                                |
C        |   By  |                                |
C        |       |                                |
C        |       |                                |
C        |       |                                |
C        Ez      |                                |
C        |       |                                |
C        |       |                                |
C        |       |                                |
C        |       *---------------Ey---------------*-----> y
C        |      /                                /
C        |     /                                /
C        |    /                                /
C        |   Ex              Bz               Ex
C        |  /                                /
C        | /                                /
C        |/                                /
C        *---------------Ey---------------*
C       /
C      /
C     /
C    x
C
C
C     Ex(i,j,k) = value of Ex at x=i+.5, y=j,    z=k
C     Ey(i,j,k) = value of Ey at x=i,    y=j+.5, z=k
C     Ez(i,j,k) = value of Ez at x=i,    y=j,    z=k+.5
C
C     Bx(i,j,k) = value of Bx at x=i,    y=j+.5, z=k+.5
C     By(i,j,k) = value of By at x=i+.5, y=j,    z=k+.5
C     Bz(i,j,k) = value of Bz at x=i+.5, y=j+.5, z=k
C
C     Maxwell's laws are implemented (to central difference
C     accuracy) in the form:
C     
C     Change of flux of B through a cell face
C         = - circulation of E around that face
C
C     Change of flux of E through a cell face of the offset grid
C         = circulation of B around that face - current through it
C
C Second half-advance of magnetic field:
      do 8 i=1,mx-1
      do 8 j=1,my-1
      do 8 k=1,mz-1
      bx(i,j,k)=bx(i,j,k) + (.5*c) *
     &          (ey(i,j,k+1)-ey(i,j,k)-ez(i,j+1,k)+ez(i,j,k))
      by(i,j,k)=by(i,j,k) + (.5*c) *
     &          (ez(i+1,j,k)-ez(i,j,k)-ex(i,j,k+1)+ex(i,j,k))
  8   bz(i,j,k)=bz(i,j,k) + (.5*c) *
     &          (ex(i,j+1,k)-ex(i,j,k)-ey(i+1,j,k)+ey(i,j,k))
C Front,right and top layers of B must be obtained from a special
C boundary routine based on Lindman's method:
      call surface(by,bz,bx,ey,ez,ex,iy,iz,ix,my,mz,mx,1,c)
      call surface(bz,bx,by,ez,ex,ey,iz,ix,iy,mz,mx,my,1,c)
      call surface(bx,by,bz,ex,ey,ez,ix,iy,iz,mx,my,mz,1,c)
C followed by post-processing the edges:
      call postedge(by,bz,bx,ey,ez,ex,iy,iz,ix,my,mz,mx,1,c)
      call postedge(bz,bx,by,ez,ex,ey,iz,ix,iy,mz,mx,my,1,c)
      call postedge(bx,by,bz,ex,ey,ez,ix,iy,iz,mx,my,mz,1,c)
C Update E-fields similarly, beginning with preliminary edge updates:
      call preledge(ey,ez,ex,by,bz,bx,-iy,-iz,-ix,my,mz,mx,lot,c)
      call preledge(ez,ex,ey,bz,bx,by,-iz,-ix,-iy,mz,mx,my,lot,c)
      call preledge(ex,ey,ez,bx,by,bz,-ix,-iy,-iz,mx,my,mz,lot,c)
C Full advance of the electric field:
      do 9 i=2,mx
      do 9 j=2,my
      do 9 k=2,mz
      ex(i,j,k)=ex(i,j,k) + c *
     &          (by(i,j,k-1)-by(i,j,k)-bz(i,j-1,k)+bz(i,j,k))
      ey(i,j,k)=ey(i,j,k) + c *
     &          (bz(i-1,j,k)-bz(i,j,k)-bx(i,j,k-1)+bx(i,j,k))
  9   ez(i,j,k)=ez(i,j,k) + c *
     &          (bx(i,j-1,k)-bx(i,j,k)-by(i-1,j,k)+by(i,j,k))
C Boundary values of the E - field must be provided at rear, left
C and bottom faces of the field domain:
      call surface(ey,ez,ex,by,bz,bx,-iy,-iz,-ix,my,mz,mx,lot,c)
      call surface(ez,ex,ey,bz,bx,by,-iz,-ix,-iy,mz,mx,my,lot,c)
      call surface(ex,ey,ez,bx,by,bz,-ix,-iy,-iz,mx,my,mz,lot,c)
C followed by post-processing of edges:
      call postedge(ey,ez,ex,by,bz,bx,-iy,-iz,-ix,my,mz,mx,lot,c)
      call postedge(ez,ex,ey,bz,bx,by,-iz,-ix,-iy,mz,mx,my,lot,c)
      call postedge(ex,ey,ez,bx,by,bz,-ix,-iy,-iz,mx,my,mz,lot,c)
C The currents due to the movement of each charge q are applied to the
C E-M fields as decrements of E-flux through cell faces. The movement
C of particles which themselves cross cell boundaries has to be split
C into several separate moves, each only within one cell. Each of
C these moves contributes to flux across twelve faces.
C Ions and electrons are processed in two loops, changing their
C attributes in-between. These loops cannot be vectorised:
C particles get processed one by one. Here is a good place to
C insert the statements for applying boundary conditions to
C the particles, such as reflection, periodicity, replacement
C by inward moving thermal or streaming particles, etc. For instance,
C to reflect particles at x=3, x=mx-2, y=3, y=my-2, z=3, z=mz-2:
C -----------
CC    loop n from 1 to ions and from maxhlf+1 to maxhlf+lecs 
CC    x0=x(n)-u(n)
CC    y0=y(n)-v(n)
CC    z0=z(n)-w(n)
CC    u(n)=u(n)*sign(1.,x(n)-3.)*sign(1.,mx-2.-x(n))
CC    x(n)=mx-2.- abs(mx-5.-abs(x(n)-3.))
CC    v(n)=v(n)*sign(1.,y(n)-3.)*sign(1.,my-2.-y(n))
CC    y(n)=my-2.- abs(my-5.-abs(y(n)-3.))
CC    w(n)=w(n)*sign(1.,z(n)-3.)*sign(1.,mz-2.-z(n))
CC    z(n)=mz-2.- abs(mz-5.-abs(z(n)-3.))
CC .. call xsplit(x(n),y(n),z(n),x0,y0,z0)
C -----------
C  Alternatively, apply periodicity to particles:
C -----------
CC    loop n from 1 to ions and from maxhlf+1 to maxhlf+lecs 
CC    x(n)=x(n)+sign(.5*(mx-5.),mx-2.-x(n))-sign(.5*(mx-5.),x(n)-3.)
CC    y(n)=y(n)+sign(.5*(my-5.),my-2.-y(n))-sign(.5*(my-5.),y(n)-3.)
CC    z(n)=z(n)+sign(.5*(mz-5.),mz-2.-z(n))-sign(.5*(mz-5.),z(n)-3.)
CC .. call xsplit(x(n),y(n),z(n),x(n)-u(n),y(n)-v(n),z(n)-w(n))
C -----------
C  The split routines call the deposit routine. They also return the verdict
C  whether the particles have left the computational domain, whether they are  
C  " in " or " .not.in ". In the solar wind application we eliminate particles
C  which have left the domain (which are  .not.in). Such particles leave their
C  charges behind on the boundary surfaces.
      q=qi
      n=1
  52  call xsplit(x(n),y(n),z(n),x(n)-u(n),y(n)-v(n),z(n)-w(n)
     &,mx,my,mz,in)
      if(in) go to 58
C  Replace by the ion from the top of the stack:
      x(n)=x(ions)
      y(n)=y(ions)
      z(n)=z(ions)
      u(n)=u(ions)
      v(n)=v(ions)
      w(n)=w(ions)
      ions=ions-1
      n=n-1
  58  n=n+1
      if(n.le.ions)go to 52
      q=qe
      n=maxhlf+1
  53  call xsplit(x(n),y(n),z(n),x(n)-u(n),y(n)-v(n),z(n)-w(n)
     &,mx,my,mz,in)
      if(in) go to 59
C  Replace by the electron from the top of the stack:
      x(n)=x(maxhlf+lecs)
      y(n)=y(maxhlf+lecs)
      z(n)=z(maxhlf+lecs)
      u(n)=u(maxhlf+lecs)
      v(n)=v(maxhlf+lecs)
      w(n)=w(maxhlf+lecs)
      lecs=lecs-1
      n=n-1
  59  n=n+1
      if(n.le.maxhlf+lecs)go to 53
C  We spread out any surface charges, and let surface electrons recombine
C  with surface ions by making the surfaces slightly conducting:
      call conduct(ix,iy,iz,mx,my,mz,ey,ez,.25)
      call conduct(iy,iz,ix,my,mz,mx,ez,ex,.25)
      call conduct(iz,ix,iy,mz,mx,my,ex,ey,.25)
C  Check if there is room for injecting new particles across the surfaces:
      if(max0(ions+212,lecs+308).gt.maxhlf)go to 68
C  ---------------
C  Inject solar wind, 180 electron-ion pairs from the front.
C  (180 per step = area*flux = area*density*vdrft = 30*30*.8*.25)
      do 61 y0=3.,my-4.,2.
      do 61 z0=3.,mz-4.5,2.5
      ions=ions+1
      x(ions)=3.+.25*rand()
      y(ions)=y0+2.*rand()
      z(ions)=z0+2.5*rand()
      u(ions)=vdrft+vth2i*(vth2i+rand()+rand()+rand()-1.5)
      v(ions)=vth2i*(rand()+rand()+rand()-1.5)
      w(ions)=vth2i*(rand()+rand()+rand()-1.5)
      lecs=lecs+1
      x(maxhlf+lecs)=x(ions)
      y(maxhlf+lecs)=y(ions)
      z(maxhlf+lecs)=z(ions)
      u(maxhlf+lecs)=vdrft+vth2e*(vth2e+rand()+rand()+rand()-1.5)
      v(maxhlf+lecs)=vth2e*(rand()+rand()+rand()-1.5)
  61  w(maxhlf+lecs)=vth2e*(rand()+rand()+rand()-1.5)
CC Inject lateral solar fluxes (32 electrons, 8 ions each face,
CC note the spacing, "spc"). Here is a case of introducing electrons and
CC ions unpaired. The code implies that for every unpaired electron an
CC immobile positive charge is created at the place of birth and stays
CC there. Likewise with each unpaired ion one creates an immobile
CC negative charge. On the lateral boundaries we create more electrons
CC than ions so these boundaries will charge up. However, since ions
CC are also less likely to escape, at least along magnetic field
CC lines or in the absence of magnetic fields, such boundary charging
CC can be expected to get neutralised by escapees. The numbers 32 and 8
CC are round-ups of 28.8 and 7.2 per face per step, accounted for as
CC area*density*vthermal/sqrt(2*pi) = 60*30*.8*.05(or .0125)/2.5
      vth2=vth2i
      spc=15.
      n=ions+1
  63  do 62 x0 = 3.,mx-3.,spc
      do 67 y0 = 3.,my-3.,spc
      x(n)=x0+spc*rand()
      y(n)=y0+spc*rand()
      z(n)=3.
      u(n)=vdrft+vth2*(rand()+rand()+rand()-1.5)
      v(n)=vth2*(rand()+rand()+rand()-1.5)
      w(n)=vth2*table(1+ifix(64.*rand()))
      n=n+1
      x(n)=x0+spc*rand()
      y(n)=y0+spc*rand()
      z(n)=mz-2.00001
      u(n)=vdrft+vth2*(rand()+rand()+rand()-1.5)
      v(n)=vth2*(rand()+rand()+rand()-1.5)
      w(n)=-vth2*table(1+ifix(64.*rand()))
  67  n=n+1
      do 62 z0 = 3.,mz-3.,spc
      x(n)=x0+spc*rand()
      z(n)=z0+spc*rand()
      y(n)=3.
      u(n)=vdrft+vth2*(rand()+rand()+rand()-1.5)
      w(n)=vth2*(rand()+rand()+rand()-1.5)
      v(n)=vth2*table(1+ifix(64.*rand()))
      n=n+1
      x(n)=x0+spc*rand()
      z(n)=z0+spc*rand()
      y(n)=my-2.00001
      u(n)=vdrft+vth2*(rand()+rand()+rand()-1.5)
      w(n)=vth2*(rand()+rand()+rand()-1.5)
      v(n)=-vth2*table(1+ifix(64.*rand()))
  62  n=n+1
      if(spc.eq.7.5)go to 64
      ions=n-1
      vth2=vth2e
      spc=7.5
      n=maxhlf+lecs+1
      go to 63
  64  lecs=n-maxhlf-1
C  Record at every timestep:
      write(*,'(3i9)')lap,ions,lecs
C ------------
C  Add the ring current to create Earth's dipole, with smoothing.
C  All current sources are smoothed. The comment "cdir$ ivdep" is
C  specifically for CRAY compilers, to encourage vectorisation.
cdir$ ivdep
  60  do 66 n=1,27
      ex(ms(n)+ie,je,ke)=ex(ms(n)+ie,je,ke)-(sm(n)*o)
      ex(ms(n)+ie,je+1,ke)=ex(ms(n)+ie,je+1,ke)+(sm(n)*o)
      ey(ms(n)+ie,je,ke)=ey(ms(n)+ie,je,ke)+(sm(n)*o)
  66  ey(ms(n)+ie+1,je,ke)=ey(ms(n)+ie+1,je,ke)-(sm(n)*o)
      if(o1.eq.0.)go to 65
C  Before the ring current 'o' has levelled off (o1=0.), we make it
C  rise gently from zero, following a cubic:
      o2=o2+o3
      o1=o1+o2
      o=o+o1
C Countdown:
  65  if(lap.ne.last)go to 6
      go to 69
  68  write(*,'(9h OVERFLOW)')
C  Record fields, particles and miscellaneous constants, for 
C  diagnostics and for continuation runs:
  69  write(7)ex,ey,ez,bx,by,bz,sm,ms,c,mx,my,mz,ix,iy,iz,lot,
     &o,o1,o2,o3,ie,je,ke,qe,qi,qme,qmi,vdrft,vth2e,vth2i,table,luck,lap
      write(8)ions,lecs,maxptl,maxhlf,x,y,z,u,v,w
      close(7)
      close(8)
      stop
      end
C -----------
      subroutine surface(bx,by,bz,ex,ey,ez,ix,iy,iz,mx,my,mz,m00,c)
      dimension bx(1),by(1),bz(1),ex(1),ey(1),ez(1)
      rs=2.*c/(1.+c)
      s=.4142136
      os=.5*(1.-s)*rs
      do 3 m=m00+iz*(mz-1),m00+iz*(mz-1)+iy*(my-2),iy
      do 1 n=m,m+ix*(mx-2),ix
   1  bz(n)=bz(n)+.5*c*(ex(n+iy)-ex(n)-ey(n+ix)+ey(n))
cdir$ ivdep
      do 3 n=m+ix,m+ix*(mx-2),ix
   3  bx(n)=bx(n)+rs*(bx(n-iz)-bx(n)+s*(bz(n)-bz(n-ix)))-os*(
     &ez(n+iy)-ez(n))-(os-c)*(ez(n+iy-iz)-ez(n-iz))-c*(ey(n)-ey(n-iz))
      do 2 m=m00+iz*(mz-1),m00+iz*(mz-1)+ix*(mx-2),ix
cdir$ ivdep
      do 4 n=m+iy,m+iy*(my-2),iy
   4  by(n)=by(n)+rs*(by(n-iz)-by(n)+s*(bz(n)-bz(n-iy)))+os*(
     &ez(n+ix)-ez(n))+(os-c)*(ez(n+ix-iz)-ez(n-iz))+c*(ex(n)-ex(n-iz))
      do 2 n=m,m+iy*(my-2),iy
   2  bz(n)=bz(n)+.5*c*(ex(n+iy)-ex(n)-ey(n+ix)+ey(n))
      return
      end
C-------------
      subroutine preledge(bx,by,bz,ex,ey,ez,ix,iy,iz,mx,my,mz,m,c)
      dimension bx(1),by(1),bz(1),ex(1),ey(1),ez(1)
      s=.4142136
      t=c/(2.*c+1.+s)
cdir$ ivdep
      do 1 n=m+iy*(my-1)+iz*(mz-1)+ix,m+iy*(my-1)+iz*(mz-1)+ix*(mx-2),ix
   1  bx(n)=bx(n-iy-iz)+(1.-4.*t)*bx(n)+(1.-2.*t)*(bx(n-iy)+bx(n-iz))
     &+s*t*(by(n)-by(n-ix)+by(n-iz)-by(n-ix-iz)
     &+bz(n)-bz(n-ix)+bz(n-iy)-bz(n-ix-iy))
      r=4./(2.*c+2.+s)
cdir$ ivdep
      do 2 n=m+iz*(mz-1),m+iz*(mz-1)+iy*(my-2),iy
   2  bx(n)=(1.-c*r)*(bx(n)+bx(n+ix))+bx(n-iz)+bx(n+ix-iz)-r*(bz(n)+c*
     &((1.-s)*(ex(n+iy)-ex(n))+(1.+s)*.25*(ez(n+iy)-ez(n)+ez(n+ix+iy)
     &-ez(n+ix)+ez(n+iy-iz)-ez(n-iz)+ez(n+ix+iy-iz)-ez(n+ix-iz))))
cdir$ ivdep
      do 3 n=m+iy*(my-1),m+iy*(my-1)+iz*(mz-2),iz
   3  bx(n)=(1.-c*r)*(bx(n)+bx(n+ix))+bx(n-iy)+bx(n+ix-iy)-r*(by(n)-c*
     &((1.-s)*(ex(n+iz)-ex(n))+(1.+s)*.25*(ey(n+iz)-ey(n)+ey(n+ix+iz)
     &-ey(n+ix)+ey(n+iz-iy)-ey(n-iy)+ey(n+ix+iz-iy)-ey(n+ix-iy))))
      p=(1.+c)*2./(1.+2.*c*(1.+c*s))
      q=c*s*2./(1.+2.*c*(1.+c*s))
cdir$ ivdep
      do 4 n=m+iy*(my-1)+iz*(mz-1),m+ix*(mx-2)+iy*(my-1)+iz*(mz-1),ix
      temp=bz(n)-.5*c*(1.-s)*(ey(n+ix)-ey(n)+ey(n+ix-iy)-ey(n-iy))
      bz(n)=bz(n-iy)-bz(n)+p*temp+q*by(n)
   4  by(n)=by(n-iz)-by(n)+p*by(n)+q*temp
      return
      end
C-------------
      subroutine postedge(bx,by,bz,ex,ey,ez,ix,iy,iz,mx,my,mz,m,c)
      dimension bx(1),by(1),bz(1),ex(1),ey(1),ez(1)
      s=.4142136
      p=(1.+c)*2./(1.+2.*c*(1.+c*s))
      q=c*s*2./(1.+2.*c*(1.+c*s))
cdir$ ivdep
      do 4 n=m+iy*(my-1)+iz*(mz-1),m+ix*(mx-2)+iy*(my-1)+iz*(mz-1),ix
      temp=by(n-iz)-.5*c*(1.-s)*(ez(n+ix)-ez(n)+ez(n+ix-iz)-ez(n-iz))
      bz(n)=bz(n-iy)+bz(n)-q*temp-p*bz(n-iy)
   4  by(n)=by(n-iz)+by(n)-q*bz(n-iy)-p*temp
      t=c/(2.*c+1.+s)
cdir$ ivdep
      do 1 n=m+iy*(my-1)+iz*(mz-1)+ix,m+iy*(my-1)+iz*(mz-1)+ix*(mx-2),ix
   1  bx(n)=bx(n)-(1.-4.*t)*bx(n-iy-iz)-(1.-2.*t)*(bx(n-iy)+bx(n-iz))
     &+s*t*(by(n)-by(n-ix)+by(n-iz)-by(n-ix-iz)
     &+bz(n)-bz(n-ix)+bz(n-iy)-bz(n-ix-iy))
      r=4./(2.*c+2.+s)
cdir$ ivdep
      do 2 n=m+iz*(mz-1),m+iz*(mz-1)+iy*(my-2),iy
   2  bx(n)=bx(n)-bx(n+ix)-(1.-c*r)*(bx(n-iz)+bx(n+ix-iz))+r*bz(n)
cdir$ ivdep
      do 3 n=m+iy*(my-1),m+iy*(my-1)+iz*(mz-2),iz
   3  bx(n)=bx(n)-bx(n+ix)-(1.-c*r)*(bx(n-iy)+bx(n+ix-iy))+r*by(n)
      return
      end
C--------
      subroutine mover(n1,n2,qm)
      common /partls/ x(110000),y(110000),z(110000),
     &u(110000),v(110000),w(110000)
C (Field components are treated as single-indexed in this subroutine)
      common /fields/ ex(79625),ey(79625),ez(79625),
     &bx(79625),by(79625),bz(79625),c,ix,iy,iz
      do 1 n=n1,n2
C Cell index & displacement in cell:
      i=x(n)
      dx=x(n)-i
      j=y(n)
      dy=y(n)-j
      k=z(n)
      dz=z(n)-k
      l=i+iy*(j-1)+iz*(k-1)
C Field interpolations are tri-linear (linear in x times linear in y
C times linear in z). This amounts to the 3-D generalisation of "area
C weighting". A modification of the simple linear interpolation formula
C         f(i+dx) = f(i) + dx * (f(i+1)-f(i))
C is needed since fields are recorded at half-integer locations in certain
C dimensions: see comments and illustration with the Maxwell part of this
C code. One then has first to interpolate from "midpoints" to "gridpoints"
C by averaging neighbors. Then one proceeds with normal interpolation.
C Combining these two steps leads to:
C   f at location i+dx  = half of f(i)+f(i-1) + dx*(f(i+1)-f(i-1))
C where now f(i) means f at location i+1/2. The halving is absorbed
C in the final scaling.
C   E-component interpolations:
      f=ex(l)+ex(l-ix)+dx*(ex(l+ix)-ex(l-ix))
      f=f+dy*(ex(l+iy)+ex(l-ix+iy)+dx*(ex(l+ix+iy)-ex(l-ix+iy))-f)
      g=ex(l+iz)+ex(l-ix+iz)+dx*(ex(l+ix+iz)-ex(l-ix+iz))
      g=g+dy*
     & (ex(l+iy+iz)+ex(l-ix+iy+iz)+dx*(ex(l+ix+iy+iz)-ex(l-ix+iy+iz))-g)
      ex0=(f+dz*(g-f))*(.25*qm)
C   ------------
      f=ey(l)+ey(l-iy)+dy*(ey(l+iy)-ey(l-iy))
      f=f+dz*(ey(l+iz)+ey(l-iy+iz)+dy*(ey(l+iy+iz)-ey(l-iy+iz))-f)
      g=ey(l+ix)+ey(l-iy+ix)+dy*(ey(l+iy+ix)-ey(l-iy+ix))
      g=g+dz*
     & (ey(l+iz+ix)+ey(l-iy+iz+ix)+dy*(ey(l+iy+iz+ix)-ey(l-iy+iz+ix))-g)
      ey0=(f+dx*(g-f))*(.25*qm)
C   ------------
      f=ez(l)+ez(l-iz)+dz*(ez(l+iz)-ez(l-iz))
      f=f+dx*(ez(l+ix)+ez(l-iz+ix)+dz*(ez(l+iz+ix)-ez(l-iz+ix))-f)
      g=ez(l+iy)+ez(l-iz+iy)+dz*(ez(l+iz+iy)-ez(l-iz+iy))
      g=g+dx*
     & (ez(l+ix+iy)+ez(l-iz+ix+iy)+dz*(ez(l+iz+ix+iy)-ez(l-iz+ix+iy))-g)
      ez0=(f+dy*(g-f))*(.25*qm)
C   ---------
C   B-component interpolations:
      f=bx(l-iy)+bx(l-iy-iz)+dz*(bx(l-iy+iz)-bx(l-iy-iz))
      f=bx(l)+bx(l-iz)+dz*(bx(l+iz)-bx(l-iz))+f+dy*
     & (bx(l+iy)+bx(l+iy-iz)+dz*(bx(l+iy+iz)-bx(l+iy-iz))-f)
      g=bx(l+ix-iy)+bx(l+ix-iy-iz)+dz*(bx(l+ix-iy+iz)-bx(l+ix-iy-iz))
      g=bx(l+ix)+bx(l+ix-iz)+dz*(bx(l+ix+iz)-bx(l+ix-iz))+g+dy*
     & (bx(l+ix+iy)+bx(l+ix+iy-iz)+dz*(bx(l+ix+iy+iz)-bx(l+ix+iy-iz))-g)
      bx0=(f+dx*(g-f))*(.125*qm/c)
C   ------------
      f=by(l-iz)+by(l-iz-ix)+dx*(by(l-iz+ix)-by(l-iz-ix))
      f=by(l)+by(l-ix)+dx*(by(l+ix)-by(l-ix))+f+dz*
     & (by(l+iz)+by(l+iz-ix)+dx*(by(l+iz+ix)-by(l+iz-ix))-f)
      g=by(l+iy-iz)+by(l+iy-iz-ix)+dx*(by(l+iy-iz+ix)-by(l+iy-iz-ix))
      g=by(l+iy)+by(l+iy-ix)+dx*(by(l+iy+ix)-by(l+iy-ix))+g+dz*
     & (by(l+iy+iz)+by(l+iy+iz-ix)+dx*(by(l+iy+iz+ix)-by(l+iy+iz-ix))-g)
      by0=(f+dy*(g-f))*(.125*qm/c)
C   ------------
      f=bz(l-ix)+bz(l-ix-iy)+dy*(bz(l-ix+iy)-bz(l-ix-iy))
      f=bz(l)+bz(l-iy)+dy*(bz(l+iy)-bz(l-iy))+f+dx*
     & (bz(l+ix)+bz(l+ix-iy)+dy*(bz(l+ix+iy)-bz(l+ix-iy))-f)
      g=bz(l+iz-ix)+bz(l+iz-ix-iy)+dy*(bz(l+iz-ix+iy)-bz(l+iz-ix-iy))
      g=bz(l+iz)+bz(l+iz-iy)+dy*(bz(l+iz+iy)-bz(l+iz-iy))+g+dx*
     & (bz(l+iz+ix)+bz(l+iz+ix-iy)+dy*(bz(l+iz+ix+iy)-bz(l+iz+ix-iy))-g)
      bz0=(f+dz*(g-f))*(.125*qm/c)
C   ---------
C   First half electric acceleration, with relativity's gamma:
      g=c/sqrt(c**2-u(n)**2-v(n)**2-w(n)**2)
      u0=g*u(n)+ex0
      v0=g*v(n)+ey0
      w0=g*w(n)+ez0
C   First half magnetic rotation, with relativity's gamma:
      g=c/sqrt(c**2+u0**2+v0**2+w0**2)
      bx0=g*bx0
      by0=g*by0
      bz0=g*bz0
      f=2./(1.+bx0*bx0+by0*by0+bz0*bz0)
      u1=(u0+v0*bz0-w0*by0)*f
      v1=(v0+w0*bx0-u0*bz0)*f
      w1=(w0+u0*by0-v0*bx0)*f
C   Second half mag. rot'n & el. acc'n:
      u0=u0+v1*bz0-w1*by0+ex0
      v0=v0+w1*bx0-u1*bz0+ey0
      w0=w0+u1*by0-v1*bx0+ez0
C   Relativity's gamma:
      g=c/sqrt(c**2+u0**2+v0**2+w0**2)
      u(n)=g*u0
      v(n)=g*v0
      w(n)=g*w0
C   Position advance:
      x(n)=x(n)+u(n)
      y(n)=y(n)+v(n)
  1   z(n)=z(n)+w(n)
      return
      end
C ------------
      subroutine conduct(ix,iy,iz,mx,my,mz,ey,ez,s)
      dimension ey(1),ez(1)
C  field arrays are treated as one-dimensional arrays in this subroutine
C  Current density components jy=s*ey, jz=s*ez are spread laterally in the
C  proportion .25,.5,.25 required by our smoothing convention for all field
C  sources.
      do 1 i=1+2*ix,1+ix*(mx-3),ix*(mx-5)
      do 1 j=i+2*iy,i+iy*(my-4),iy
      do 1 k=j+2*iz,j+iz*(mz-4),iz
      temp=.25*s*ey(k)
      ey(k-ix)=ey(k-ix)-temp
      ey(k+ix)=ey(k+ix)-temp
      ey(k)=ey(k)-2.*temp
      temp=.25*s*ez(k)
      ez(k-ix)=ez(k-ix)-temp
      ez(k+ix)=ez(k+ix)-temp
   1  ez(k)=ez(k)-2.*temp
      return
      end
C  --------
      subroutine xsplit(x,y,z,x0,y0,z0,mx,my,mz,in)
      logical in
      in=.true.
      if((ifix(x).ne.ifix(x0)).and.((x-x0).ne.0.))go to 1
      call ysplit(x,y,z,x0,y0,z0,mx,my,mz,in)
      return
  1   x1=.5*(1+ifix(x)+ifix(x0))
      y1=y0+(y-y0)*((x1-x0)/(x-x0))
      z1=z0+(z-z0)*((x1-x0)/(x-x0))
      call ysplit(x1,y1,z1,x0,y0,z0,mx,my,mz,in)
      if(.not.in)return
      in=(x1.gt.3.0).and.(x1.lt.mx-2.0)
      if(in)call ysplit(x,y,z,x1,y1,z1,mx,my,mz,in)
      return
      end
C --------
      subroutine ysplit(x,y,z,x0,y0,z0,mx,my,mz,in)
      logical in
      if((ifix(y).ne.ifix(y0)).and.((y-y0).ne.0.))go to 1
      call zsplit(x,y,z,x0,y0,z0,mx,my,mz,in)
      return
  1   y1=.5*(1+ifix(y)+ifix(y0))
      z1=z0+(z-z0)*((y1-y0)/(y-y0))
      x1=x0+(x-x0)*((y1-y0)/(y-y0))
      call zsplit(x1,y1,z1,x0,y0,z0,mx,my,mz,in)
      if(.not.in)return
      in=(y1.gt.3.0).and.(y1.lt.my-2.0)
      if(in)call zsplit(x,y,z,x1,y1,z1,mx,my,mz,in)
      return
      end
C --------
      subroutine zsplit(x,y,z,x0,y0,z0,mx,my,mz,in)
      logical in
      if((ifix(z).ne.ifix(z0)).and.((z-z0).ne.0.))go to 1
      call depsit(x,y,z,x0,y0,z0)
      return
  1   z1=.5*(1+ifix(z)+ifix(z0))
      x1=x0+(x-x0)*((z1-z0)/(z-z0))
      y1=y0+(y-y0)*((z1-z0)/(z-z0))
      call depsit(x1,y1,z1,x0,y0,z0)
      in=(z1.gt.3.0).and.(z1.lt.mz-2.0)
      if(in)call depsit(x,y,z,x1,y1,z1)
      return
      end
C ----------------------
      subroutine depsit(x,y,z,x0,y0,z0)
C (Field components are treated as single-indexed in this subroutine)
      common /fields/ ex(79625),ey(79625),ez(79625),
     &bx(79625),by(79625),bz(79625),c,ix,iy,iz
      common /spread/q,sm(27),ms(27)
C   cell indices of half-way point:
      i=.5*(x+x0)
      j=.5*(y+y0)
      k=.5*(z+z0)
C   displacements in cell of half-way point:
      dx=.5*(x+x0) - i
      cx=1.-dx
      dy=.5*(y+y0) - j
      cy=1.-dy
      dz=.5*(z+z0) - k
      cz=1.-dz
      l=i+iy*(j-1)+iz*(k-1)
C   current elements:
      qu=q*(x-x0)
      qv=q*(y-y0)
      qw=q*(z-z0)
      delt=.08333333*qu*(y-y0)*(z-z0)
C  Directive specifically for the CRAY cft77 compiler:
cdir$ ivdep
C  (This will make the compiler use the "gather-scatter" hardware.)
C ---------------
C  If one desires NO smoothing (risking the presence of alias-prone
C  high harmonics), one can replace the statement "do 1 n=1,27" by
C     n=14
C  and boost the value of  q  by a factor 8.
C ---------------
!$acc parallel loop
      do 1 n=1,27
      m=ms(n)+l
      s=sm(n)*delt
C
      su=sm(n)*qu
      ex(m+iy+iz)=ex(m+iy+iz)-su*dy*dz-s
      ex(m+iz)=ex(m+iz)-su*cy*dz+s
      ex(m+iy)=ex(m+iy)-su*dy*cz+s
      ex(m)=ex(m)-su*cy*cz-s
C
      sv=sm(n)*qv
      ey(m+iz+ix)=ey(m+iz+ix)-sv*dz*dx-s
      ey(m+ix)=ey(m+ix)-sv*cz*dx+s
      ey(m+iz)=ey(m+iz)-sv*dz*cx+s
      ey(m)=ey(m)-sv*cz*cx-s
C
      sw=sm(n)*qw
      ez(m+ix+iy)=ez(m+ix+iy)-sw*dx*dy-s
      ez(m+iy)=ez(m+iy)-sw*cx*dy+s
      ez(m+ix)=ez(m+ix)-sw*dx*cy+s
  1   ez(m)=ez(m)-sw*cx*cy-s
!$acc end parallel loop
      return
      end
C===================== end of tristan.f ==========================
