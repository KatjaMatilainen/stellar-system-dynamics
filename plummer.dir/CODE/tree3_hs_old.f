C***********************************************************************
C
C
                             PROGRAM treecode
C
C
C                      Version 3: April 1, 1990
C
C                  Lars Hernquist, Princeton University.
C
C
C***********************************************************************
C
C
C     A code to evolve self-gravitating systems using the hierarchical
C     tree method developed by Barnes and Hut (Nature 324, 446 [1986])
C     and implemented in FORTRAN by Hernquist (Ap. J. Suppl. 64, 715
C     [1987]; Comp. Phys. Comm. 48, 107 [1988]).  This version has 
C     been optimized for supercomputers and is fully vectorized.  The 
C     code is written in standard FORTRAN, although CRAY-specific 
C     vector intrinsic functions have been used (see below).
C     
C     In this version, vectorization of the tree walks is achieved by 
C     simultaneously processing all cells at the same level in the 
C     tree, as discussed by Hernquist (J. Comp. Phys, 87, 137 [1990]).  
C     The gravitational force calculation and nearest neighbor 
C     search proceed for a single particle at a time, in serial order.
C
C     The gravitational field is smoothed according to a cubic spline
C     kernel.  The kernel is computed by linear interpolation from a 
C     look-up table.
C
C     A self-starting leap-frog integrator is used, as described by
C     Hernquist and Katz (Ap. J. Suppl., 70, 419 [1989]).
C     
C     The gravitational softening length of each particle is
C     optionally variable.  Softening lengths are determined from the 
C     requirement that each softening volume contain a constant number 
C     of near neighbors.  The softening length for cells is always 
C     taken to be a mass-weighted mean of the softening lengths for 
C     all the particles it contains.
C
C     The computational system of units is determined by the input
C     data, with the assumption that G=1 .  Particles are not
C     required to have identical masses.
C
C     This file contains the complete source code for CRAYs.  In
C     order to run on other machines the file treeutil.f must also
C     be linked along with either treeutilsun.f (UNIX machines) or
C     treeutilvax.f (VMS machines).  Common-block variables are
C     defined in the include file treedefs.h.
C
C     Two input data files are required to run this code: a parameter
C     file, which is read in through the subroutine inparams, and a
C     body data file, read by subroutine inbods.  Both are ASCII and
C     their structure is defined in the subroutines which read them.
C     Three output files are created: an ASCII log file, an ASCII
C     body data file containing the final state of the system, and
C     a binary body data file.  The log and binary body files are
C     updated every noutlog and noutbod steps, respectively.
C
C     WARNINGS -- To avoid excessive overhead, noutlog should be
C                 larger than 1, typically ~ 10, depending on the
C                 number of steps.
C
C                 When compiling on VAX's, avoid using optimization.
C
C     Please report all problems or suggestions for improvements to
C     this code to lars@rigel.princeton.edu.
C
C
C=======================================================================
C
C
C     This is the top-level evolution program treecode.  Its tasks are:
C
C          1) to initialize file structures and global variables;
C          2) to input parameters and the initial system state;
C          3) to advance the state of the system for a given number
C             of timesteps;
C          4) to perform a diagnostic analysis of the system at
C             each time step (energy, angular momentum, etc.);
C          5) to periodically record the state of the system;
C          6) and to terminate the simulation and close data files.
C
C
C=======================================================================
C
C
C     Basic global variables/parameters:
C
C          acc         : acceleration components of a body.
C          acsmooth    : table of smoothed gravitational acceleration.
C          cellsize    : linear sizes of cells.
C          cputime     : cpu time (secs) used during the simulation.
C          cputime0    : cumulative cpu time at start of run.
C          cputime1    : cumulative cpu time at end of run.
C          dgravsum    : option to compute gravity by direct 
C                        summation (.TRUE.).
C          dtime       : the timestep.
C          dtime2      : timestep/2.
C          eps         : gravitational smoothing parameter.
C          epsvect     : values of gravitational softening length
C                        for each particle and each cell.
C          ektot       : total system kinetic energy.
C          eptot       : total system gravitational potential energy.
C          esofttot    : energy correction to virial theorem from
C                        particle softening.
C          etot        : total energy of the system.
C          four        : the constant 4.
C          headline    : identification string for the run.
C          incells     : number of cells currently in use.
C          incellsg    : number of cells used to compute gravity.
C          inpteps     : option to read in values of gravitational
C                        softening parameter from particle data
C                        file -- for collisionless particles only.
C          intfour     : the integer constant 4.
C          intone      : the integer constant 1.
C          intthree    : the integer constant 3.
C          inttwo      : the integer constant 2.
C          intzero     : the integer constant 0.
C          log2        : the constant log10(2).
C          mass        : masses of bodies and cells.
C          minusone    : the constant -1.
C          minustwo    : the constant -2.
C          mtot        : total mass of the system.
C          nbodies     : total number of bodies.
C          nbodsmax    : maximum number of bodies.
C          ncells      : maximum number of cells.
C          ndim        : number of spatial dimensions.
C          ninterp     : number of values in look-up tables.
C          noutbod     : frequency of system record outputs.
C          noutlog     : frequency of outputs to log file.
C          nsavg       : average particle number in softening volume.
C          nsmax       : maximum particle number in softening volume.
C          nsmin       : minimum particle number in softening volume.
C          nstot       : average particle number in softening volume.
C          nsteps      : total number of timesteps in the simulation.
C          nsubcell    : number of subcells per cell.
C          nsvolume    : number of collisionless particles per
C                        softening volume; used only if variabls
C                        is .TRUE.
C          nsvtol      : fractional tolerance in number of neighbors
C                        relative to nsvolume.  A real number 
C                        typically ~ 0.05.
C          ntavg       : average length of interaction lists.
C          ntmax       : largest interaction list in current time step.
C          ntmin       : shortest interaction list in current time step.
C          nttot       : sum of interaction lists in current time step.
C          one         : the constant 1.
C          onehalf     : the constant 0.5.
C          phi         : gravitational potential.
C          phsmooth    : table of look-up values for grav. potential.
C          pos         : coordinates of bodies, c.m. coords. of cells.
C          quad        : quadrupole moments of cells.
C          restart     : option to restart run from a SYSDUMP file.
C          rmin        : coords. of lower-left corner of system box.
C          root        : pointer to the top of the tree.
C          rsize       : length of the system box.
C          selfgrav    : option to turn off (.FALSE.) system self-
C                        gravity.
C          subp        : pointers to descendents of a cell.
C          three       : the constant 3.
C          tiny        : a small number used to prevent divergences.
C          tnow        : current system time.
C          tol         : accuracy parameter.
C          tol2        : tol * tol.
C          tpos        : current position time.
C          ttree       : time of last update of gravitational tree.
C          two         : the constant 2.
C          usequad     : option to use (.TRUE.) quadrupole terms.
C          variabls    : option to use (.TRUE.) variable gravitational
C                        softening lengths for collisionless particles.
C          vel         : velocity components of a body.
C          zero        : the constant 0.
C          zero02      : the constant 0.02.
C
C-----------------------------------------------------------------------
C
C   Definitions specific to input/output.
C
C          ireclog                        : log file record counter.
C          uterm, upars, ulog, ubodsin,   : logical i/o unit numbers.
C            ubodsout,uboddump,utermfil,
C            ucrash,uindump
C          parsfile, logfile, ibodfile,   : character names of files.
C            obodfile,dumpfile,termfile,
C            crashfil,indumpf
C
C-----------------------------------------------------------------------
C
C   Definitions specific to individual particle time steps.
C
C          endstep     : indicator that a full step is complete.
C          etol        : tolerance parameter to select new time step.
C          inittbin    : initial time step bin of all particles.          
C          itimestp    : defines time step of each particle by 
C                        dtime/itimestp.
C          mintstep    : defines minimum allowed timestep by 
C                        dtime/mintstep.
C          npactive    : number of particles needing acceleration.
C          ntvector    : minimum number of particles in each time
C                        step bin.
C          otimestp    : previous value of itimestp.
C          pactive     : indices of particles requiring acceleration.
C          stime       : fraction of large timestep remaining.
C          tsteppos    : time step with which to step positions.
C          upbin       : last time step bin which was updated.
C
C-----------------------------------------------------------------------
C
C   Definitions specific to vectorized tree construction, vectorized 
C   tree walk, and vectorized tree search for nearest neighbors.
C
C          asubp       : subpointers for active bodies or cells.
C          bodlist     : list of active bodies (i.e. not yet leaves).
C          bottom      : coordinates of bottom edges of cells.
C          celllist    : list of cells.
C          groupbod    : list of bodies in search groups.
C          groups      : list of grouped cells.
C          ingroup     : number of particles in groups for searching.
C          isubset     : indices of subsets of active bodies or cells.
C          ngroups     : number of cells containing ingroup particles.
C          npercell    : number of particles in a cell.
C          nnavg       : average number of near neighbors.
C          nnmax       : maximum number of near neighbors.
C          nnmin       : minimum number of near neighbors.
C          nnear       : number of neighbors.
C          nnearlis    : number of neighbors in neighbor lists.
C          nntot       : total number of near neighbors.
C          nworkvec    : length of temporary work array workvect.  It
C                        should be set to 9*max length of grav
C                        interaction list.
C          parent      : parents of active bodies or cells.
C          pgroupb     : pointer to list of bodies in each search group.
C          pnear       : pointer to list of nearest neighbors.
C          subindex    : subindices for active bodies or cells.
C          subpvect    : vector equivalenced to subp.
C          templist    : temporary vector to swap arrays.
C          tempvect    : temporary storage vector.
C          workvect    : temporary work array.
C
C
C=======================================================================
C
C
C     Data structure used to compute gravitational field:
C
C          The body/cell tree structure is assumed to be of the
C          form discussed by Barnes and Hut.  Schematically, for
C          three dimensions (i.e. eight subcells per cell):
C
C         +-------------------------------------------------+
C  root-->| CELL:  mass, pos, quad, /, o, /, /, /, /, o, /  |
C         +----------------------------|--------------|-----+
C                                      |              |
C     +--------------------------------+              |
C     |                                               |
C     |   +----------------------------------+        |
C     +-->| BODY:  mass, pos, vel, acc, phi  |        |
C         +----------------------------------+        |
C                                                     |
C     +-----------------------------------------------+
C     |
C     |   +--------------------------------------------------+
C     +-->| CELL:  mass, pos, quad,  o, /, /, o, /, /, o, /  |
C         +--------------------------|--------|--------|-----+
C                                    |        |        |
C                                   etc.     etc.     etc.
C
C
C          The body/cell information is stored in arrays which
C          incorporate both bodies and cells.  For physical
C          quantities relevant to both bodies and cells, such as
C          mass and position, the array indices range from
C          1 --> nbodsmax + ncells.  For those quantities defined
C          only for bodies, such as velocity, the array indices
C          range from 1 --> nbodsmax.  For information relevant
C          to cells alone, such as pointers to descendants, the
C          array indices range from nbodsmax + 1 --> nbodsmax +
C          ncells.  With this convention, pointers can refer to
C          either bodies or cells without conflict.
C
C          The identification of a given unit as a body or a cell
C          is determined by the pointer to the body/cell.  For a
C          body, p is less than or equal to nbodsmax, while for a
C          cell, p > nbodsmax.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER n

C=======================================================================

C   Initialize state of the system.
C   -------------------------------
        CALL initsys
C            -------

C   Advance system state for a given number of steps.
C   -------------------------------------------------

        DO 100 n=1,nsteps

           CALL stepsys(n)
C               -------

 100    CONTINUE

C   Terminate the simulation.
C   -------------------------
        CALL endrun
C            ------

        STOP
        END
C***********************************************************************
C
C
                        SUBROUTINE accgrav(option)
C
C
C***********************************************************************
C
C
C     Subroutine to compute the gravitational acceleration for all of
C     the bodies.  Vectorization is achieved by processing all of the
C     cells at a given level in the tree simultaneously.  The local
C     variable option indicates whether the code is to compute the
C     potential and/or acceleration.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        CHARACTER*4 option
        INTEGER p,i,j,nterms

C=======================================================================

C   Initialize the interaction list diagnostics.
C   --------------------------------------------
        nttot=0
        ntmin=nbodies
        ntmax=0

C   Main loop over all bodies.
C   --------------------------

        DO 100 i=1,npactive

           p=pactive(i)

C   Establish interaction lists.
C   ----------------------------
           IF(.NOT.dgravsum) THEN

              CALL treewalk(p,nterms)
C                  --------
           ELSE

              nterms=nbodies

              DO 50 j=1,nbodies
                 bodlist(j)=j
 50           CONTINUE

           ENDIF

C   Compute potential and acceleration.
C   -----------------------------------

           CALL gravsum(p,nterms,option)
C               -------
C   Update diagnostics, subtracting self-interaction term.
C   ------------------------------------------------------
           nterms=nterms-1
           nttot=nttot+nterms
           IF(nterms.LT.ntmin) ntmin=nterms
           IF(nterms.GT.ntmax) ntmax=nterms

 100    CONTINUE

C   Compute average number of force terms per body.
C   -----------------------------------------------
        ntavg=nttot/npactive

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE celledge
C
C
C***********************************************************************
C
C
C     Subroutine to initialize coordinates of edges of cells for 
C     nearest neighbor searching.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER isgnsub(8,3),i,j,k,nclist,nclist2,nsubc,indlist

        SAVE isgnsub

        DATA isgnsub/4*0,4*1,0,0,1,1,0,0,1,1,0,1,0,1,0,1,0,1/

C=======================================================================

C   Initialize properties of root cell.
C   -----------------------------------

        DO 10 k=1,ndim
           bottom(root,k)=rmin(k)
 10     CONTINUE

C   Place root on list of cells.
C   ----------------------------
        celllist(1)=root
        nclist=1

C   Loop until all cells are processed.
C   -----------------------------------

 200    CONTINUE

        IF(nclist.GT.0) THEN

C   Select non-zero subcells from celllist.
C   ---------------------------------------
           nclist2=0

           DO 40 j=1,nsubcell

              DO 20 i=1,nclist
                 asubp(i)=subp(celllist(i),j)
 20           CONTINUE

              CALL WHENIGT(nclist,asubp,1,nbodsmax,isubset,nsubc)
  
              nclist2=nclist2+nsubc

              IF(nclist2.GT.ncells.OR.nclist2.GT.nbodsmax)
     &           CALL terror(' array overflow in celledge ')
C                     ------

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 30 i=1,nsubc
                 indlist=nclist2-nsubc+i
                 templist(indlist)=asubp(isubset(i))
                 parent(indlist)=celllist(isubset(i))
                 subindex(indlist)=j
 30           CONTINUE

 40        CONTINUE

           nclist=nclist2
        
CVD$ NODEPCHK
CDIR$ IVDEP
           DO 60 i=1,nclist
              celllist(i)=templist(i)
              bottom(celllist(i),1)=bottom(parent(i),1)+
     &                      isgnsub(subindex(i),1)*cellsize(celllist(i))    
              bottom(celllist(i),2)=bottom(parent(i),2)+
     &                      isgnsub(subindex(i),2)*cellsize(celllist(i))    
              bottom(celllist(i),3)=bottom(parent(i),3)+
     &                      isgnsub(subindex(i),3)*cellsize(celllist(i))    
 60        CONTINUE

           GO TO 200
        
        ENDIF 

        DO 70 i=1,nbodies
           bottom(i,1)=pos(i,1)
           bottom(i,2)=pos(i,2)
           bottom(i,3)=pos(i,3)
 70     CONTINUE

        RETURN
        END
C***********************************************************************
C
C
                           SUBROUTINE celleps
C
C
C***********************************************************************
C
C
C     Subroutine to compute gravitational softening lengths of
C     cells, processing cells  in order of increasing size.  The 
C     permutation vector is stored in the common variable celllist.
C     Vectorization is achieved by simultaneously processing all 
C     cells at the same level in the hierarchy.  The softening length 
C     for each cell is a mass-weighted mean of the softening lengths 
C     of its descendents.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p,fcell,lcell,i,j,nnodes,lcf1

C=======================================================================
        
C   Generate permutation of cells, according to cellsize.
C   -----------------------------------------------------

        DO 5 i=1,incells
           celllist(i)=nbodsmax+incells-(i-1)
 5      CONTINUE

C   Initialize properties of cells.
C   -------------------------------

        DO 10 p=nbodsmax+1,nbodsmax+incells
           epsvect(p)=0.0
 10     CONTINUE

C   Process cells in order of increasing size.
C   ------------------------------------------

        fcell=1

 40     CONTINUE

        IF(fcell.LE.incells) THEN

C   Determine which cells to process.
C   ---------------------------------

           DO 50 i=fcell,incells
              IF(ABS(cellsize(celllist(i))-cellsize(celllist(fcell)))
     &           .LT.0.01*cellsize(celllist(fcell))) THEN

                 lcell=i
              ELSE
                 GO TO 60
              ENDIF
 50        CONTINUE                    

 60        CONTINUE

           lcf1=lcell-fcell+1

           IF(lcf1.GT.ncells.OR.lcf1.GT.nbodsmax)
     &        CALL terror(' lcf1 overflow in celleps ')
C                  ------

C   Compute properties of the selected cells, looping over subcells.
C   ----------------------------------------------------------------

           DO 110 j=1,nsubcell

              DO 70 i=fcell,lcell
                 asubp(i-fcell+1)=subp(celllist(i),j)
 70           CONTINUE

              CALL WHENIGT(lcell-fcell+1,asubp,1,0,isubset,nnodes)

              IF(nnodes.GT.ncells.OR.nnodes.GT.nbodsmax)
     &           CALL terror(' array overflow in celleps ')
C                     ------

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 80 i=1,nnodes
                 parent(i)=celllist(isubset(i)+fcell-1)
                 asubp(i)=subp(parent(i),j)
                 epsvect(parent(i))=epsvect(parent(i))+mass(asubp(i))*
     &                              epsvect(asubp(i))
 80           CONTINUE

 110       CONTINUE


CVD$ NODEPCHK
CDIR$ IVDEP
           DO 120 i=fcell,lcell
              epsvect(celllist(i))=epsvect(celllist(i))/
     &           mass(celllist(i))
 120       CONTINUE

           fcell=lcell+1

           GO TO 40

        ENDIF

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE cellnumb
C
C
C***********************************************************************
C
C
C     Subroutine to compute the number of particles per cell for the
C     nearest neighbor search.  Vectorization is achieved by 
C     simultaneously processing all cells at the same level in the
C     hierarchy.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p,fcell,lcell,i,j,nnodes,nsubc,nsubb

C=======================================================================

C   Generate permutation of cells, according to cellsize.
C   -----------------------------------------------------

        DO 5 i=1,incells
           celllist(i)=nbodsmax+incells-(i-1)
 5      CONTINUE

C   Initialize cell properties.
C   ---------------------------

        DO 10 p=nbodsmax+1,nbodsmax+incells
           npercell(p)=0
 10     CONTINUE

C   Process cells in order of increasing size.
C   ------------------------------------------

        fcell=1

 20     CONTINUE

        IF(fcell.LE.incells) THEN

C   Determine which cells to process.
C   ---------------------------------

           DO 30 i=fcell,incells
              IF(ABS(cellsize(celllist(i))-cellsize(celllist(fcell)))
     &           .LT.0.01*cellsize(celllist(fcell))) THEN

                 lcell=i
              ELSE
                 GO TO 40
              ENDIF
 30        CONTINUE

 40        CONTINUE

C   Compute properties of selected cells, looping over subcells.
C   ------------------------------------------------------------

           DO 100 j=1,nsubcell

              DO 50 i=fcell,lcell
                 asubp(i-fcell+1)=subp(celllist(i),j)
 50           CONTINUE

              CALL WHENIGT(lcell-fcell+1,asubp,1,0,isubset,nnodes)

              IF(nnodes.GT.ncells.OR.nnodes.GT.nbodsmax)
     &           CALL terror(' array overflow in cellnumb ')
C                     ------

              DO 60 i=1,nnodes
                 parent(i)=celllist(isubset(i)+fcell-1)
                 asubp(i)=subp(parent(i),j)
 60           CONTINUE

              CALL WHENIGT(nnodes,asubp,1,nbodsmax,isubset,nsubc)               

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 70 i=1,nsubc
                 templist(i)=parent(isubset(i))
                 npercell(templist(i))=npercell(templist(i))+
     &              npercell(asubp(isubset(i)))
 70           CONTINUE

              CALL WHENILE(nnodes,asubp,1,nbodsmax,isubset,nsubb)

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 80 i=1,nsubb
                 templist(i)=parent(isubset(i))
                 npercell(templist(i))=npercell(templist(i))+1
 80           CONTINUE

 100       CONTINUE

           fcell=lcell+1

           GO TO 20

        ENDIF

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE checkinp
C
C
C***********************************************************************
C
C
C     Subroutine to check consistency of input parameters and data,
C     output warnings to the terminal and/or log file, and terminate
C     the simulation if necessary.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C=======================================================================

        IF(nsteps.LT.0.OR.nsteps.GT.10000)
     &     CALL terror(' input error for parameter nsteps ')
C               ------

        IF(noutbod.LT.0)
     &     CALL terror(' input error for parameter noutbod ')
C               ------

        IF(noutlog.LT.0)
     &     CALL terror(' input error for parameter noutlog ')
C               ------

        IF(dtime.LE.0.0.OR.dtime.GT.1.e20)
     &     CALL terror(' input error for parameter dtime ')
C               ------

        IF((.NOT.dgravsum).AND.(tol.LT.0.0.OR.tol.GT.1.5))
     &     CALL terror(' input error for parameter tol ')
C               ------

        IF(eps.LT.0.0.OR.eps.GT.1.e20)
     &     CALL terror(' input error for parameter eps ')
C               ------

        IF(variabls.AND.(nsvolume.LE.0.OR.nsvolume.GT.nbodsmax))
     &     CALL terror(' input error for parameter nsvolume ')
C               ------

        IF(inittbin.LT.1.OR.inittbin.GT.262144.OR.inittbin.GT.mintstep)
     &     CALL terror(' input error for parameter inittbin ')
C               ------

        IF(mintstep.LT.1.OR.mintstep.GT.262144)
     &     CALL terror(' input error for parameter mintstep ')
C               ------

        IF(etol.LE.0.0.OR.etol.GT.1.e20)
     &     CALL terror(' input error for parameter etol ')
C               ------

        IF(ntvector.LE.0.OR.ntvector.GT.262144)
     &     CALL terror(' input error for parameter ntvector ')
C               ------

        RETURN
        END
C***********************************************************************
C
C
                         SUBROUTINE corrpos(rc)
C
C
C***********************************************************************
C
C
C     Subroutine to apply a correction factor to the positions to
C     maintain second order accuracy when outputting particle data
C     to body data file or when computing energy diagnostics.  The
C     argument rc indicates whether the correction factor is to be
C     applied (correct) or removed (reset).
C  
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        CHARACTER*7 rc
        INTEGER p,k
        REAL dt2,rcsign

C=======================================================================

C   Loop over all spatial coordinates for all bodies.
C   -------------------------------------------------
  
        IF(rc.EQ.'correct') THEN
           rcsign=-1.
        ELSE
           rcsign=1.
        ENDIF

        DO 200 k=1,ndim
           DO 100 p=1,nbodies
              dt2=(dtime/itimestp(p))**2
              pos(p,k)=pos(p,k)+rcsign*acc(p,k)*dt2/8.
 100       CONTINUE
 200    CONTINUE

        RETURN
        END

C***********************************************************************
C
C
                           SUBROUTINE endrun
C
C
C***********************************************************************
C
C
C     Subroutine to end the simulation.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        REAL second

C=======================================================================

        cputime1=SECOND()

        CALL outcpu
C            ------
        CALL stopout
C            -------

        RETURN
        END
C***********************************************************************
C
C
                            SUBROUTINE energy
C
C
C***********************************************************************
C
C
C     Subroutine to compute diagnostics for the system: total energy,
C     total kinetic energy, total potential energy, angular momentum,
C     center of mass coordinates, and center of mass velocity.  The
C     local variable p is a pointer to the bodies.  The local
C     variables cmpos, cmvel, and amvec are center of mass position 
C     and velocity, and total angular momentum of the system, 
C     respectively.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p,k
        REAL cmpos(ndim),cmvel(ndim),amvec(3)

C=======================================================================

C   Zero the accumulators for system diagnostics.
C   ---------------------------------------------
        mtot=0.
        ektot=0.
        eptot=0.

        DO 100 k=1,ndim
           cmpos(k)=0.
           cmvel(k)=0.
 100    CONTINUE

        DO 120 k=1,3
           amvec(k)=0.
 120    CONTINUE

C-----------------------------------------------------------------------
C   Loop over bodies to compute system mass and potential energy.
C-----------------------------------------------------------------------

        DO 150 p=1,nbodies
           mtot=mtot+mass(p)
           eptot=eptot+.5*mass(p)*phi(p)
 150    CONTINUE

C-----------------------------------------------------------------------
C   Compute system kinetic energy, components of center of mass
C   position and velocity.
C-----------------------------------------------------------------------

        DO 250 k=1,ndim
           DO 200 p=1,nbodies
              ektot=ektot+.5*mass(p)*vel(p,k)*vel(p,k)
              cmpos(k)=cmpos(k)+mass(p)*pos(p,k)
              cmvel(k)=cmvel(k)+mass(p)*vel(p,k)
 200       CONTINUE
           cmvel(k)=cmvel(k)/mtot
           cmpos(k)=cmpos(k)/mtot
 250    CONTINUE

C   Compute total system energy.
C   ----------------------------
        etot=ektot+eptot

C-----------------------------------------------------------------------
C   Compute angular momentum of the system.
C-----------------------------------------------------------------------

        IF(ndim.EQ.2) THEN

           DO 300 p=1,nbodies
              amvec(3)=amvec(3)+mass(p)*(pos(p,1)*vel(p,2)-
     &                 pos(p,2)*vel(p,1))
 300       CONTINUE

        ELSE IF(ndim.EQ.3) THEN

           DO 400 p=1,nbodies
              amvec(1)=amvec(1)+mass(p)*(pos(p,2)*vel(p,3)-
     &                 pos(p,3)*vel(p,2))
              amvec(2)=amvec(2)+mass(p)*(pos(p,3)*vel(p,1)-
     &                 pos(p,1)*vel(p,3))
              amvec(3)=amvec(3)+mass(p)*(pos(p,1)*vel(p,2)-
     &                 pos(p,2)*vel(p,1))
 400       CONTINUE

        ENDIF

C   Write diagnostics to the log file.
C   ----------------------------------
        CALL outenrgy(amvec,cmpos,cmvel)
C            --------

        RETURN
        END
C***********************************************************************
C
C
                  SUBROUTINE enlargee(pbody,nnearbod)
C
C
C***********************************************************************
C
C
C     Subroutine to increase softening length of body pbody so that
C     the softening volume contains roughly nsvolume neighbors.  The
C     neighbor list is returned from the subroutine nearpeps in the
C     vector nearlist, which is equivalenced to the common array
C     parent.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER i,pbody,nhtest,nearlist(nbodsmax),nnearbod,
     &          jsubset(nbodsmax),npnear,keepnear(nbodsmax)
        REAL hmin,hmax,htest,fourhsm2,testnn,hsearch,dx,dy,dz

        EQUIVALENCE (nearlist(1),parent(1)),(jsubset(1),subindex(1)),
     &              (keepnear(1),asubp(1))

C=======================================================================

        hsearch=2.0
        npnear=0

 5      CONTINUE

        IF(npnear.LT.nsvolume) THEN

           hsearch=1.5*hsearch

           IF(hsearch.GT.1.e5)
     &        CALL terror(' search error in enlargee ')
C                  ------

           CALL nearpeps(hsearch,pbody,npnear)
C               --------

           GO TO 5

        ENDIF
        
        nhtest=npnear

        DO 10 i=1,npnear
           dx=pos(pbody,1)-pos(nearlist(i),1)
           dy=pos(pbody,2)-pos(nearlist(i),2)
           dz=pos(pbody,3)-pos(nearlist(i),3)
           tempvect(i)=dx**2+dy**2+dz**2
           isubset(i)=i
 10     CONTINUE

        testnn=ABS(REAL(npnear-nsvolume))/REAL(nsvolume)

        IF(npnear.LE.nsvolume.OR.testnn.LE.nsvtol) THEN

           htest=hsearch*epsvect(pbody)

        ELSE

           hmin=0.
           hmax=hsearch*epsvect(pbody)

 20        CONTINUE

           IF(ABS(REAL(nhtest-nsvolume)).GT.nsvtol*REAL(nsvolume)) THEN

              htest=0.5*(hmin+hmax)
              fourhsm2=four*htest*htest

              CALL WHENFLT(npnear,tempvect,1,fourhsm2,isubset,nhtest)

              IF(nhtest.GT.nsvolume) THEN
                 hmax=htest
              ELSE
                 hmin=htest
              ENDIF

              GO TO 20

           ENDIF

        ENDIF

        nstot=nstot+nhtest
        IF(nhtest.LT.nsmin) nsmin=nhtest
        IF(nhtest.GT.nsmax) nsmax=nhtest
        nnear(pbody)=nhtest
        epsvect(pbody)=htest

        RETURN
        END
C***********************************************************************
C
C
                      SUBROUTINE findnear(p,npnear)
C
C
C***********************************************************************
C
C
C     Subroutine to search for nearest neighbors of the grouped
C     cell p.  Vectorization is achieved by processing all cells
C     at a given level in the tree simultaneously.  The list of
C     neighbors is returned through the vector nearlist, which is
C     equivalenced to the common array bodlist.
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p,i,nsubdiv,nearlist(nbodsmax),npnear,nnodes,nkeep,
     &          nodelist(ncells),keepnear(nbodsmax),keepstak(nbodsmax)
        LOGICAL testcrit
        REAL pbottom(ndim),ptop(ndim),dx,dy,dz,xnode,ynode,znode,
     &       groupx,groupy,groupz

        EQUIVALENCE (nearlist(1),bodlist(1)),(nodelist(1),celllist(1)),
     &              (keepnear(1),parent(1)),(keepstak(1),asubp(1))

C=======================================================================

        npnear=0

        CALL srchbox(p,npnear,pbottom,ptop)
C            -------

        nnodes=intone
        nodelist(1)=root

 20     CONTINUE

        IF(nnodes.GT.0) THEN

CVD$ NODEPCHK
CDIR$ IVDEP
           DO 30 i=1,nnodes
              xnode=bottom(nodelist(i),1)+0.5*cellsize(nodelist(i))
              ynode=bottom(nodelist(i),2)+0.5*cellsize(nodelist(i))
              znode=bottom(nodelist(i),3)+0.5*cellsize(nodelist(i))
              groupx=0.5*(ptop(1)+pbottom(1))
              groupy=0.5*(ptop(2)+pbottom(2))
              groupz=0.5*(ptop(3)+pbottom(3))
              dx=0.0
              dy=0.0
              dz=0.0
              testcrit=     pbottom(1).LE.(bottom(nodelist(i),1)+
     &                                   cellsize(nodelist(i))+dx)
     &                 .AND.pbottom(2).LE.(bottom(nodelist(i),2)+
     &                                   cellsize(nodelist(i))+dy)
     &                 .AND.pbottom(3).LE.(bottom(nodelist(i),3)+
     &                                   cellsize(nodelist(i))+dz)
     &                 .AND.ptop(1).GE.(bottom(nodelist(i),1)+dx)
     &                 .AND.ptop(2).GE.(bottom(nodelist(i),2)+dy)
     &                 .AND.ptop(3).GE.(bottom(nodelist(i),3)+dz)
     &                 .AND.nodelist(i).NE.groups(p)
              IF(testcrit.AND.nodelist(i).LE.nbodsmax) THEN
                 keepnear(i)=2
              ELSE
                 keepnear(i)= -2
              ENDIF
              IF(testcrit.AND.nodelist(i).GT.nbodsmax) THEN
                 keepstak(i)=2
              ELSE
                 keepstak(i)= -2
              ENDIF
 30        CONTINUE

           CALL WHENIGT(nnodes,keepnear,1,0,isubset,nkeep)

           IF(npnear+nkeep.GT.nbodsmax)
     &        CALL terror(' array overflow in findnear ')
C                  ------

           DO 40 i=1,nkeep
              nearlist(npnear+i)=nodelist(isubset(i))
 40        CONTINUE

           npnear=npnear+nkeep

           CALL WHENIGT(nnodes,keepstak,1,0,isubset,nsubdiv)

           IF(8*nsubdiv.GT.nbodsmax.OR.8*nsubdiv.GT.ncells)
     &        CALL terror(' asubp overflow in findnear ')
C                  ------

CVD$ NODEPCHK
CDIR$ IVDEP
           DO 50 i=1,nsubdiv
              asubp(i)=subp(nodelist(isubset(i)),1)
              asubp(i+nsubdiv)=subp(nodelist(isubset(i)),2)
              asubp(i+2*nsubdiv)=subp(nodelist(isubset(i)),3)
              asubp(i+3*nsubdiv)=subp(nodelist(isubset(i)),4)
              asubp(i+4*nsubdiv)=subp(nodelist(isubset(i)),5)
              asubp(i+5*nsubdiv)=subp(nodelist(isubset(i)),6)
              asubp(i+6*nsubdiv)=subp(nodelist(isubset(i)),7)
              asubp(i+7*nsubdiv)=subp(nodelist(isubset(i)),8)
 50        CONTINUE

           CALL WHENNE(8*nsubdiv,asubp,1,0,isubset,nnodes)

           DO 60 i=1,nnodes
              nodelist(i)=asubp(isubset(i))
 60        CONTINUE

           GO TO 20

        ENDIF

        RETURN
        END
C***********************************************************************
C
C
                       SUBROUTINE gravity(option)
C
C
C***********************************************************************
C
C
C     Subroutine to compute gravitational potential and acceleration.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        CHARACTER*4 option

C=======================================================================

        IF(selfgrav) THEN

           IF((.NOT.dgravsum).OR.variabls) CALL maketree
C                                               --------
           IF(variabls) THEN

              CALL maketeps
C                  --------
              CALL stepeps
C                  -------
           ENDIF

           IF(.NOT.dgravsum) CALL celleps
C                                 -------
           CALL accgrav(option)
C               -------

        ENDIF

        RETURN
        END
C***********************************************************************
C
C
                  SUBROUTINE gravsum(p,nterms,option)
C
C
C***********************************************************************
C
C
C     Subroutine to compute the monopole and quadrupole contributions
C     to the potential and acceleration components for body p.  The
C     interaction list is contained in the vector iterms, which is
C     equivalenced to the common array bodlist.
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER maxnterm

        PARAMETER(maxnterm=nworkvec/9)

        CHARACTER*4 option
        INTEGER p,i,qindex(nbodsmax),qterms(nbodsmax),smindex(nbodsmax),
     &          nterms,iterms(nbodsmax),nqterms
        REAL r3inveff(maxnterm),rinveff(maxnterm),drdeldrg,pmass,
     &       drdotdr(maxnterm),phsm,drsm,accsm,dx(maxnterm),
     &       dy(maxnterm),dz(maxnterm),qr5inv(maxnterm),acci,
     &       phiquad(maxnterm),sdrdotdr,r2inveff(maxnterm),esoftsum

        EQUIVALENCE (iterms(1),bodlist(1)),(qindex(1),templist(1)),
     &              (qterms(1),isubset(1)),(smindex(1),parent(1)),
     &              (dx(1),workvect(1)),(dy(1),workvect(maxnterm+1)),
     &              (dz(1),workvect(2*maxnterm+1)),(r3inveff(1),
     &              workvect(3*maxnterm+1)),(rinveff(1),
     &              workvect(4*maxnterm+1)),(drdotdr(1),
     &              workvect(5*maxnterm+1)),(qr5inv(1),
     &              workvect(6*maxnterm+1)),(phiquad(1),
     &              workvect(7*maxnterm+1)),(r2inveff(1),
     &              workvect(8*maxnterm+1))

C=======================================================================

        IF(nterms.GT.maxnterm)
     &     CALL terror(' array overflow in gravsum ')
C               ------

C-----------------------------------------------------------------------
C   Compute monopole contribution; temporarily set mass of body p to
C   zero to avoid possible self-interaction contribution.
C-----------------------------------------------------------------------

        pmass=mass(p)
        mass(p)=0.

C   Loop over interaction list.
C   ---------------------------

CVD$ NODEPCHK
CDIR$ IVDEP
        DO 30 i=1,nterms
           dx(i)=pos(p,1)-pos(iterms(i),1)
           dy(i)=pos(p,2)-pos(iterms(i),2)
           dz(i)=pos(p,3)-pos(iterms(i),3)
           drdotdr(i)=dx(i)**2+dy(i)**2+dz(i)**2+tiny*
     &                (epsvect(p)+epsvect(iterms(i)))**2/4.
           sdrdotdr=SQRT(drdotdr(i))
           rinveff(i)=1./sdrdotdr
           r3inveff(i)=rinveff(i)/drdotdr(i)
           drdeldrg=sdrdotdr*ninterp/(epsvect(p)+epsvect(iterms(i)))
           smindex(i)=drdeldrg
           IF(ninterp.LT.smindex(i)) smindex(i)=ninterp
           IF(one.LT.drdeldrg-smindex(i)) THEN
              drsm=one
           ELSE
              drsm=drdeldrg-smindex(i)
           ENDIF
           phsm=(1.-drsm)*phsmooth(smindex(i))+
     &          drsm*phsmooth(1+smindex(i))
           accsm=(1.-drsm)*acsmooth(smindex(i))+
     &           drsm*acsmooth(1+smindex(i))
           rinveff(i)=phsm*rinveff(i)
           r3inveff(i)=accsm*r3inveff(i)
 30     CONTINUE

        IF(option.NE.'acc ') THEN

           IF(npactive.NE.nbodies) THEN

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 40 i=1,nterms
                 phi(p)=phi(p)-mass(iterms(i))*rinveff(i)
 40           CONTINUE

           ELSE

              esoftsum=0.0

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 45 i=1,nterms
                 phi(p)=phi(p)-mass(iterms(i))*rinveff(i)
                 esoftsum=esoftsum+mass(iterms(i))*(rinveff(i)-
     &                    drdotdr(i)*r3inveff(i))
 45           CONTINUE

           ENDIF

        ENDIF

        IF(option.NE.'pot ') THEN

CVD$ NODEPCHK
CDIR$ IVDEP
           DO 50 i=1,nterms
              acci=mass(iterms(i))*r3inveff(i)
              acc(p,1)=acc(p,1)-dx(i)*acci
              acc(p,2)=acc(p,2)-dy(i)*acci
              acc(p,3)=acc(p,3)-dz(i)*acci
 50        CONTINUE

        ENDIF

C   Reset mass of body p.
C   ---------------------
        mass(p)=pmass
 
        IF(npactive.EQ.nbodies) esofttot=esofttot+0.5*mass(p)*esoftsum

C   If required, compute quadrupole contribution.
C   ---------------------------------------------
        IF(usequad) THEN

C   Filter out bodies.
C   ------------------

           CALL WHENIGT(nterms,iterms,1,nbodsmax,qindex,nqterms)
 
C   Compute quadrupole interaction from cells.
C   ------------------------------------------
CVD$ NODEPCHK
CDIR$ IVDEP
           DO 60 i=1,nqterms
              qterms(i)=iterms(qindex(i))
              r2inveff(i)=rinveff(qindex(i))*rinveff(qindex(i))
              qr5inv(i)=r3inveff(qindex(i))*r2inveff(i)
              phiquad(i)=(-.5*((dx(qindex(i))**2-dz(qindex(i))**2)*
     &              quad(qterms(i),1)+(dy(qindex(i))**2-
     &              dz(qindex(i))**2)*quad(qterms(i),4))-
     &              (dx(qindex(i))*dy(qindex(i))*quad(qterms(i),2)+
     &              dx(qindex(i))*dz(qindex(i))*quad(qterms(i),3)+
     &              dy(qindex(i))*dz(qindex(i))*quad(qterms(i),5)))*
     &              qr5inv(i)
 60        CONTINUE

           IF(option.NE.'acc ') THEN

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 70 i=1,nqterms
                 phi(p)=phi(p)+phiquad(i)
 70           CONTINUE

           ENDIF

           IF(option.NE.'pot ') THEN

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 80 i=1,nqterms
                 phiquad(i)=5.*phiquad(i)*r2inveff(i)
                 acc(p,1)=acc(p,1)+dx(qindex(i))*phiquad(i)+
     &                  (dx(qindex(i))*quad(qterms(i),1)+
     &                  dy(qindex(i))*quad(qterms(i),2)+
     &                  dz(qindex(i))*quad(qterms(i),3))*qr5inv(i)
                 acc(p,2)=acc(p,2)+dy(qindex(i))*phiquad(i)+
     &                  (dy(qindex(i))*quad(qterms(i),4)+
     &                  dx(qindex(i))*quad(qterms(i),2)+
     &                  dz(qindex(i))*quad(qterms(i),5))*qr5inv(i)
                 acc(p,3)=acc(p,3)+dz(qindex(i))*phiquad(i)+
     &                  (dz(qindex(i))*(-quad(qterms(i),1)-
     &                  quad(qterms(i),4))+dx(qindex(i))*
     &                  quad(qterms(i),3)+dy(qindex(i))*
     &                  quad(qterms(i),5))*qr5inv(i)
 80           CONTINUE

           ENDIF

        ENDIF
 
        RETURN
        END
C***********************************************************************
C
C
                         SUBROUTINE groupcel
C
C
C***********************************************************************
C
C
C     Subroutine to locate cells containing ingroup particles, subject
C     to the constraint that the smoothing lengths of the particles
C     do not span too large a range.  Lists of such cells and the 
C     particles that lie within them are created. 
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER pstack,stack(nbodsmax),spointer,i,j,iblist,gpointer,
     &          gstack(nbodsmax),pgstack,numgroup,ismin,ismax
        LOGICAL testsubd
        REAL hmax,hmin

        EQUIVALENCE (stack(1),parent(1)),(gstack(1),asubp(1))

C=======================================================================

C   Locate cells containing ingroup particles.
C   ------------------------------------------
        stack(1)=root
        spointer=1
        ngroups=0
        iblist=0

 120    CONTINUE

        pstack=stack(spointer)
        spointer=spointer-1

C   If ingroup exceeded, subdivide cell; otherwise record cell.
C   -----------------------------------------------------------

        IF(pstack.LE.nbodsmax) THEN
           ngroups=ngroups+1
           groups(ngroups)=pstack
           iblist=iblist+1
           groupbod(iblist)=pstack
           pgroupb(ngroups)=iblist
        ELSE
           IF(npercell(pstack).GT.ingroup) THEN
              DO 130 j=1,nsubcell
                 IF(subp(pstack,j).GT.0) THEN
                    spointer=spointer+1
                    IF(spointer.GT.nbodsmax)
     &                 CALL terror(' array overflow in groupcel ')
C                           ------
                    stack(spointer)=subp(pstack,j)
                 ENDIF
 130          CONTINUE
           ELSE

C   Walk through tree, listing bodies within grouped cells.
C   -------------------------------------------------------

              gpointer=1
              gstack(1)=pstack
              numgroup=0

 150          CONTINUE

              pgstack=gstack(gpointer)
              gpointer=gpointer-1

              IF(pgstack.GT.nbodsmax) THEN

                 DO 160 j=1,nsubcell
                    IF(subp(pgstack,j).GT.0) THEN
                       gpointer=gpointer+1
                       IF(gpointer.GT.nbodsmax)
     &                    CALL terror(' array overflow in groupcel ')
C                              ------
                       gstack(gpointer)=subp(pgstack,j)
                    ENDIF
 160             CONTINUE

              ELSE

                 numgroup=numgroup+1
                 bodlist(numgroup)=pgstack

              ENDIF

              IF(gpointer.GT.0) GO TO 150

C   Subdivide if range in smoothing or softening lengths too large.
C   ---------------------------------------------------------------

              DO 170 i=1,numgroup
                 tempvect(i)=epsvect(bodlist(i))
 170          CONTINUE

              hmax=tempvect(ISMAX(numgroup,tempvect,1))
              hmin=tempvect(ISMIN(numgroup,tempvect,1))

              testsubd=cellsize(pstack).GT.5.*hmax.AND.hmax.GT.0.
              testsubd=testsubd.OR.hmax.GT.2.*hmin

              IF(testsubd) THEN

                 DO 180 j=1,nsubcell
                    IF(subp(pstack,j).GT.0) THEN
                       spointer=spointer+1
                       IF(spointer.GT.nbodsmax)
     &                    CALL terror(' array overflow in groupcel ')
C                              ------
                       stack(spointer)=subp(pstack,j)
                    ENDIF
 180             CONTINUE

              ELSE

                 ngroups=ngroups+1
                 groups(ngroups)=pstack
                 pgroupb(ngroups)=iblist+1

                 DO 190 i=1,numgroup
                    iblist=iblist+1
                    groupbod(iblist)=bodlist(i)
 190             CONTINUE

              ENDIF

           ENDIF
        ENDIF

        IF(spointer.GT.0) GO TO 120

        pgroupb(ngroups+1)=iblist+1

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE hackcell
C
C
C***********************************************************************
C
C
C     Subroutine to compute masses, center of mass coordinates,
C     and optional quadrupole moments of cells, processing cells
C     in order of increasing size.  The permutation vector is
C     stored in the common variable celllist.  Vectorization is
C     achieved by simultaneously processing all cells at the
C     same level in the hierarchy.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p,fcell,lcell,i,j,k,l,m,n,nsubb,nsubc,nnodes,mupper,
     &          lcf1

C=======================================================================
        
C   Generate permutation of cells, according to cellsize.
C   -----------------------------------------------------

        DO 5 i=1,incells
           celllist(i)=nbodsmax+incells-(i-1)
 5      CONTINUE

C   Initialize properties of cells.
C   -------------------------------

        DO 10 p=nbodsmax+1,nbodsmax+incells
           mass(p)=0.
           npercell(p)=0
           pos(p,1)=0.
           pos(p,2)=0.
           pos(p,3)=0.
 10     CONTINUE

        IF(usequad) THEN
           DO 30 k=1,2*ndim-1
              DO 20 p=nbodsmax+1,nbodsmax+incells
                 quad(p,k)=0.
 20           CONTINUE
 30        CONTINUE
        ENDIF

C   Process cells in order of increasing size.
C   ------------------------------------------

        fcell=1

 40     CONTINUE

        IF(fcell.LE.incells) THEN

C   Determine which cells to process.
C   ---------------------------------

           DO 50 i=fcell,incells
              IF(ABS(cellsize(celllist(i))-cellsize(celllist(fcell)))
     &           .LT.0.01*cellsize(celllist(fcell))) THEN

                 lcell=i
              ELSE
                 GO TO 60
              ENDIF
 50        CONTINUE                    

 60        CONTINUE

           lcf1=lcell-fcell+1

           IF(lcf1.GT.ncells.OR.lcf1.GT.nbodsmax)
     &        CALL terror(' lcf1 overflow in hackcell ')
C                  ------

C   Compute properties of the selected cells, looping over subcells.
C   ----------------------------------------------------------------

           DO 110 j=1,nsubcell

              DO 70 i=fcell,lcell
                 asubp(i-fcell+1)=subp(celllist(i),j)
 70           CONTINUE

              CALL WHENIGT(lcell-fcell+1,asubp,1,0,isubset,nnodes)

              IF(nnodes.GT.ncells.OR.nnodes.GT.nbodsmax)
     &           CALL terror(' array overflow in hackcell ')
C                     ------

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 80 i=1,nnodes
                 parent(i)=celllist(isubset(i)+fcell-1)
                 asubp(i)=subp(parent(i),j)
                 mass(parent(i))=mass(parent(i))+mass(asubp(i))
                 pos(parent(i),1)=pos(parent(i),1)+mass(asubp(i))*
     &                            pos(asubp(i),1)
                 pos(parent(i),2)=pos(parent(i),2)+mass(asubp(i))*
     &                            pos(asubp(i),2)
                 pos(parent(i),3)=pos(parent(i),3)+mass(asubp(i))*
     &                            pos(asubp(i),3)
 80           CONTINUE

              CALL WHENIGT(nnodes,asubp,1,nbodsmax,isubset,nsubc)

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 90 i=1,nsubc
                 templist(i)=parent(isubset(i))
                 npercell(templist(i))=npercell(templist(i))+
     &              npercell(asubp(isubset(i)))
 90           CONTINUE

              CALL WHENILE(nnodes,asubp,1,nbodsmax,isubset,nsubb)

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 100 i=1,nsubb
                 templist(i)=parent(isubset(i))
                 npercell(templist(i))=npercell(templist(i))+1
 100          CONTINUE

 110       CONTINUE

CVD$ NODEPCHK
CDIR$ IVDEP
           DO 120 i=fcell,lcell
              pos(celllist(i),1)=pos(celllist(i),1)/mass(celllist(i))
              pos(celllist(i),2)=pos(celllist(i),2)/mass(celllist(i))
              pos(celllist(i),3)=pos(celllist(i),3)/mass(celllist(i))
 120       CONTINUE

C   Compute optional quadrupole moments.
C   ------------------------------------

           IF(usequad) THEN

              DO 210 j=1,nsubcell

                 DO 130 i=fcell,lcell
                    asubp(i-fcell+1)=subp(celllist(i),j)
 130             CONTINUE

                 CALL WHENIGT(lcell-fcell+1,asubp,1,0,isubset,nnodes)

CVD$ NODEPCHK
CDIR$ IVDEP
                 DO 140 i=1,nnodes
                    parent(i)=celllist(isubset(i)+fcell-1)
                    asubp(i)=subp(parent(i),j)
 140             CONTINUE

                 CALL WHENIGT(nnodes,asubp,1,nbodsmax,isubset,nsubc)

                 IF(ndim.GT.2) THEN
                    mupper=2
                 ELSE
                    mupper=ndim
                 ENDIF

                 DO 200 m=1,mupper
                    DO 190 n=m,ndim

                       l=(m-1)*(ndim-1)+n

CVD$ NODEPCHK
CDIR$ IVDEP
                       DO 150 i=1,nnodes
                          quad(parent(i),l)=quad(parent(i),l)+
     &                       mass(asubp(i))*(3.*(pos(asubp(i),m)-
     &                       pos(parent(i),m))*(pos(asubp(i),n)-
     &                       pos(parent(i),n)))
 150                   CONTINUE

                       IF(m.EQ.n) THEN
                          DO 170 k=1,ndim
CVD$ NODEPCHK
CDIR$ IVDEP
                             DO 160 i=1,nnodes
                                quad(parent(i),l)=quad(parent(i),l)-
     &                             mass(asubp(i))*(pos(asubp(i),k)-
     &                             pos(parent(i),k))**2
 160                         CONTINUE
 170                      CONTINUE
                       ENDIF

CVD$ NODEPCHK
CDIR$ IVDEP
                       DO 180 i=1,nsubc
                          templist(i)=parent(isubset(i))
                          quad(templist(i),l)=quad(templist(i),l)+
     &                       quad(asubp(isubset(i)),l)
 180                   CONTINUE

 190                CONTINUE
 200             CONTINUE

 210          CONTINUE

           ENDIF

           fcell=lcell+1
             
           GO TO 40

        ENDIF

        RETURN
        END

C***********************************************************************
C
C
                           SUBROUTINE inbods
C
C
C***********************************************************************
C
C
C     Subroutine to read in the data associated with the bodies.  The
C     records are assumed to be of the form: nbodies, ndim, time,
C     mass(1)...mass(n), x(1)...x(n), y(1)...y(n), ..., vx(1)...vx(n),
C     vy(1)...vy(n), ....
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p,ndimi

C=======================================================================
 
        OPEN(UNIT=ubodsin,FILE=ibodfile,STATUS='OLD')

C   Read in body data.
C   ------------------

        READ(ubodsin,*) nbodies
        READ(ubodsin,*) ndimi
        READ(ubodsin,*) tnow

        IF(nbodies.GT.nbodsmax.OR.ndimi.NE.ndim)
     &     CALL terror(' error in inbods--inconsistent inputs ')
C               ------

        DO 10 p=1,nbodies
           READ(ubodsin,*) mass(p)
 10     CONTINUE

        DO 20 p=1,nbodies
           READ(ubodsin,*) pos(p,1),pos(p,2),pos(p,3)
 20     CONTINUE

        DO 30 p=1,nbodies
           READ(ubodsin,*) vel(p,1),vel(p,2),vel(p,3)
 30     CONTINUE

        IF(inpteps) THEN
           DO 40 p=1,nbodies
              READ(ubodsin,*) epsvect(p)
 40        CONTINUE
        ENDIF

        CLOSE(ubodsin)
 
        RETURN
        END
C***********************************************************************
C
C
                            SUBROUTINE indump
C
C
C***********************************************************************
C
C
C     Subroutine to read in system dump from an ascii data file.
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p

C=======================================================================
 
        OPEN(UNIT=uindump,FILE=indumpf,STATUS='OLD')

C   Output system state.
C   --------------------

        READ(uindump,*) nbodies

        READ(uindump,*) tnow,tpos
        READ(uindump,*) dtime
 
        READ(uindump,*) mtot
        READ(uindump,*) ektot
        READ(uindump,*) eptot
        READ(uindump,*) esofttot

        READ(uindump,*) upbin
        READ(uindump,*) stime
        READ(uindump,*) tsteppos
        READ(uindump,*) endstep

        DO 70 p=1,nbodies
           READ(uindump,*) mass(p)
 70     CONTINUE

        DO 80 p=1,nbodies
           READ(uindump,*) pos(p,1),pos(p,2),pos(p,3)
 80     CONTINUE

        DO 90 p=1,nbodies
           READ(uindump,*) vel(p,1),vel(p,2),vel(p,3)
 90     CONTINUE

        DO 105 p=1,nbodies
           READ(uindump,*) epsvect(p)
 105    CONTINUE

        DO 210 p=1,nbodies
           READ(uindump,*) nnear(p)
 210    CONTINUE

        DO 220 p=1,nbodies
           READ(uindump,*) phi(p)
 220    CONTINUE

        DO 230 p=1,nbodies
           READ(uindump,*) acc(p,1),acc(p,2),acc(p,3)
 230    CONTINUE

        DO 240 p=1,nbodies
           READ(uindump,*) itimestp(p)
 240    CONTINUE

        DO 250 p=1,nbodies
           READ(uindump,*) otimestp(p)
 250    CONTINUE

        CLOSE(UNIT=uindump)

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE initeps
C
C
C***********************************************************************
C
C
C     Subroutine to initialize gravitational softening lengths for
C     collisionless particles, allowing for a target number of near 
C     neighbors.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER i,j,nnearfix
        REAL avgdens,epsavg

C=======================================================================

        IF(variabls.OR.eps.EQ.zero) THEN

           DO 10 i=1,nbodies
              avgdens=nbodies/rsize**3
              epsvect(i)=0.01*(nsvolume/(4.*3.141592654*avgdens/3.))**
     &                   (one/3.)
 10        CONTINUE

           DO 40 j=1,20

              CALL neighcol('predict')
C                  --------

              DO 30 i=1,nbodies
                 IF(nnear(i).EQ.0) THEN
                    nnearfix=1
                 ELSE
                    nnearfix=0
                 ENDIF
                 epsvect(i)=0.5*epsvect(i)*((REAL(nsvolume)/
     &                      REAL(nnear(i)+nnearfix))**(one/3.)+one)
 30           CONTINUE

 40        CONTINUE           

           IF(.NOT.variabls) THEN

              epsavg=0.

              DO 50 i=1,nbodies
                 epsavg=epsavg+epsvect(i)
 50           CONTINUE

              epsavg=epsavg/nbodies

              DO 60 i=1,nbodies
                 epsvect(i)=epsavg
 60           CONTINUE

           ENDIF

        ELSE

           DO 70 i=1,nbodies
              epsvect(i)=eps
 70        CONTINUE

        ENDIF           

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE initpars
C
C
C***********************************************************************
C
C
C     Subroutine to initialize system parameters that depend on
C     either the input data or defined PARAMETERS.  The local
C     variable p is a pointer to the bodies/cells.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p,i
        REAL xw,xw2,deldrg,xw3,xw4

C=======================================================================

C   Initialize misc. useful numbers.
C   --------------------------------
        minusone = -1.
        minustwo = -2.
        tiny=1.e-20
        zero=0.
        zero02=0.02
        intzero=0
        intone=1
        inttwo=2
        intthree=3
        intfour=4
        one=1.
        two=2.
        three=3.
        four=4.
        log2=LOG(two)
        onehalf=0.5

C   Initialize position and velocity times, 1/2 timestep.
C   -----------------------------------------------------
        IF(.NOT.restart) THEN
           tpos=tnow
        ENDIF

        ttree=tnow-1.
        dtime2=.5*dtime
        tol2=tol*tol

C   Initialize size parameter for bodies.
C   -------------------------------------

        DO 5 p=1,nbodies
           cellsize(p)=0.
 5      CONTINUE

        rsize=0.
        rmin(1)=0.
        rmin(2)=0.
        rmin(3)=0.

        IF(.NOT.restart) THEN

           DO 7 p=1,nbodies
              nnear(p)=0
 7         CONTINUE

        ENDIF

C-----------------------------------------------------------------------
C   Initialize variables and arrays for gravitational field smoothing 
C   interpolation.  Interpolation performed in distance.
C-----------------------------------------------------------------------
        deldrg=2./ninterp

        phsmooth(0)=7.*SQRT(tiny)/5.
        acsmooth(0)=4.*SQRT(tiny)*tiny/3.

        DO 30 i=1,1+ninterp
           xw=i*deldrg
           xw2=xw*xw
           xw3=xw2*xw
           xw4=xw2*xw2
           IF(xw.LE.one) THEN
              phsmooth(i)=-2.*xw3*(one/3.-3.*xw2/20.+xw3/20.)+7.*xw/5.
              acsmooth(i)=xw3*(4./3.-6.*xw2/5.+0.5*xw3)
           ELSE
              phsmooth(i)=-one/15.+8.*xw/5.-xw3*(4./3.-xw+3.*xw2/10.-
     &                    xw3/30.)
              acsmooth(i)=-one/15.+8.*xw3/3.-3.*xw4+6.*xw3*xw2/5.-
     &                    xw4*xw2/6.
           ENDIF
           IF(xw.GE.two) THEN
              phsmooth(i)=one
              acsmooth(i)=one
           ENDIF
 30     CONTINUE

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE initpos
C
C
C***********************************************************************
C
C
C     Subroutine to initialize the positions of the bodies for the
C     initial timestep, dtime/initbin.  The local variable p is a
C     pointer to the bodies.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p,k
        REAL dt2,dt

C=======================================================================

C   Loop over all spatial coordinates for all bodies.
C   -------------------------------------------------
  
        dt=dtime/inittbin
        dt2=dt**2

        DO 200 k=1,ndim
           DO 100 p=1,nbodies
              pos(p,k)=pos(p,k)+acc(p,k)*dt2/8.
 100       CONTINUE
 200    CONTINUE

        RETURN
        END

C***********************************************************************
C
C
                          SUBROUTINE initstep
C
C
C***********************************************************************
C
C
C     Subroutine to initialize parameters dealing with individual
C     particle time steps.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER i

C=======================================================================

        DO 5 i=1,nbodies
           itimestp(i)=inittbin
           otimestp(i)=inittbin
 5      CONTINUE

        upbin=mintstep
        stime=0.0
        endstep=.TRUE.
        npactive=nbodies
        tsteppos=0.

        DO 55 i=1,npactive
           pactive(i)=i
 55     CONTINUE

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE initsys
C
C
C***********************************************************************
C
C
C     Subroutine to initialize the state of the system.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        REAL second

C=======================================================================

C   Begin timing.
C   -------------
        cputime0=SECOND()

C-----------------------------------------------------------------------
C   Open data files, read input parameters and initial system state,
C   and initialize system parameters.
C-----------------------------------------------------------------------
        CALL startout
C            --------
        CALL inparams
C            --------

        IF(.NOT.restart) THEN

           CALL inbods
C               ------
        ELSE

           CALL indump
C               ------
        ENDIF

        CALL checkinp
C            --------
        CALL initpars
C            --------

        IF(.NOT.restart) CALL initstep
C                             --------
        IF(.NOT.restart) THEN

C   Zero out potential and acceleration.
C   ------------------------------------
           CALL zeroacc
C               -------
           CALL zeropot
C               -------

C   Initialize gravitational softening lengths.
C   -------------------------------------------
           IF(variabls.OR.eps.EQ.0.0) CALL maketeps
C                                          --------
           IF(.NOT.inpteps) CALL initeps
C                                -------
           IF(variabls) CALL neighcol('correct')
C                            --------

C   Compute gravitational potential and acceleration.
C   -------------------------------------------------
           CALL gravity('both')
C               -------

        ENDIF

C   Output system state.
C   --------------------
        IF(.NOT.restart) THEN
           CALL outstate(0)
C               --------
        ELSE
           CALL outhead
C               -------
        ENDIF

        IF(.NOT.restart) CALL initpos
C                             -------

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE inparams
C
C
C***********************************************************************
C
C
C     Subroutine to read in parameters.
C
C     Input parameters:
C
C        headline  : identification string for the run.
C        nsteps    : number of timesteps.
C        noutbod   : output system state once every nsteps/noutbod 
C                    steps.
C        noutlog   : output logfile data once every nsteps/noutlog
C                    steps.
C        dtime     : the timestep.
C        tol       : error tolerance; 0.0 => exact (PP) calculation.
C        eps       : potential softening parameter.
C        inpteps   : option to read gravitational softening lengths
C                    of collisionless particles from input data.
C        variabls  : option to use variable gravitational softening
C                    lengths for collisionless particles.
C        nsvolume  : number of collisionless particles per softening
C                    volume.
C        nsvtol    : fractional tolerance in number of neighbors
C                    relative to nsvolume; typically ~ 0.05.
C        usequad   : option to include (.TRUE.) quadrupole terms.
C        dgravsum  : option to use direct summation gravity (.TRUE.).
C        restart   : option to restart run from a SYSDUMP file.
C-----------------------------------------------------------------------
C        inittbin  : initial time step bin for all particles.
C        mintstep  : parameter defining minimum allowed time step.
C        etol      : energy tolerance to select timestep.
C        ntvector  : minumum number of particles per time step bin.
C-----------------------------------------------------------------------
C        selfgrav  : option to turn off (.FALSE.) system self-gravity.
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        CHARACTER *1 pcomment

C=======================================================================
 
        OPEN(UNIT=upars,FILE=parsfile,STATUS='OLD')

C   Read parameters, close the file.
C   --------------------------------

        READ(upars,'(a)') pcomment

        READ(upars,'(a)') headline
        READ(upars,*) nsteps
        READ(upars,*) noutbod
        READ(upars,*) noutlog
        READ(upars,*) dtime
        READ(upars,*) tol
        READ(upars,*) eps
        READ(upars,*) inpteps
        READ(upars,*) variabls
        READ(upars,*) nsvolume
        READ(upars,*) nsvtol
        READ(upars,*) usequad
        READ(upars,*) dgravsum
        READ(upars,*) restart

        READ(upars,'(a)') pcomment

        READ(upars,*) inittbin
        READ(upars,*) mintstep
        READ(upars,*) etol
        READ(upars,*) ntvector

        READ(upars,'(a)') pcomment

        READ(upars,*) selfgrav

c  hs
        read(upars,*) noutani,niraf,biraf


        CLOSE(UNIT=upars)
 
        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE loadtree
C
C
C***********************************************************************
C
C
C     Subroutine to insert the bodies into the tree.  The process is 
C     vectorized over active bodies.  Active bodies are those which 
C     are not yet in place in the tree, as leaves.  The local variables
C     pm1 and nindex are used to convert back and forth between 
C     physical coordinates and subcell coordinates.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER k,p,nindex(ndim),j,i,nbodlist,nclist,nclist2,
     &          nsubset,indcell,nsubbod1,nbodtemp,iupper,ilower,
     &          ncl2nsub
        REAL pm1(nsubcell,ndim)

        SAVE nindex,pm1

        DATA pm1/4*-1.,4*1.,2*-1.,2*1.,2*-1.,2*1.,-1.,1.,-1.,1.,
     &            -1.,1.,-1.,1./,nindex/4,2,1/

C=======================================================================

C   Deallocate old tree, compute coordinates of center of root cell.
C   ----------------------------------------------------------------
        incells=1
        root=nbodsmax+1

        DO 5 j=1,nsubcell
           subp(root,j)=0
 5      CONTINUE

        cellsize(root)=rsize

        DO 10 k=1,ndim
           pos(root,k)=rmin(k)+0.5*rsize
           bottom(root,k)=rmin(k)
 10     CONTINUE

C-----------------------------------------------------------------------
C   Place all bodies on active body list, having root as parent; place
C   root on active cell list.
C-----------------------------------------------------------------------

        ilower=1
        iupper=nbodies

        DO 20 i=ilower,iupper
           parent(i-ilower+1)=root
           bodlist(i-ilower+1)=i
 20     CONTINUE

        nbodlist=iupper-ilower+1
        celllist(1)=root
        nclist=1

C   Loop until no bodies are left active.
C   -------------------------------------

 200    CONTINUE

        IF(nclist.GT.0) THEN

C   Compute subindices for all active bodies.
C   -----------------------------------------
           DO 30 i=1,nbodlist
              subindex(i)=1
 30        CONTINUE

           DO 50 k=1,ndim
              DO 40 i=1,nbodlist
                 IF(pos(bodlist(i),k).GE.pos(parent(i),k)) 
     &                  subindex(i)=subindex(i)+nindex(k)
 40           CONTINUE
 50        CONTINUE

C   Compute number of bodies in each subcell.
C   -----------------------------------------
           DO 60 i=1,nbodlist
              subp(parent(i),subindex(i))=subp(parent(i),subindex(i))+1
 60        CONTINUE

C-----------------------------------------------------------------------
C   Open all subcells with more than one body, placing them on active 
C   cell list.
C-----------------------------------------------------------------------
           nclist2=0

           DO 110 j=1,nsubcell

              DO 70 i=1,nclist
                 asubp(i)=subp(celllist(i),j)
 70           CONTINUE

              CALL WHENIGT(nclist,asubp,1,1,isubset,nsubset)

              incells=incells+nsubset

              IF(incells.GT.ncells.OR.incells.GT.nbodsmax) 
     &           CALL terror(' overflow in loadtree')
C                     ------
              indcell=incells-nsubset+nbodsmax

              ncl2nsub=nclist2+nsubset

              IF(ncl2nsub.GT.nbodsmax.OR.ncl2nsub.GT.ncells)
     &           CALL terror(' nclist2 overflow in loadtree ')
C                     ------

              DO 90 k=1,nsubcell
                 DO 80 i=1,nsubset
                    subp(indcell+i,k)=0
 80              CONTINUE
 90           CONTINUE

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 100 i=1,nsubset
                 p=indcell+i
                 asubp(i)=celllist(isubset(i))
                 subp(asubp(i),j)=p
                 cellsize(p)=cellsize(asubp(i))*0.5
                 templist(nclist2+i)=p
                 pos(p,1)=pos(asubp(i),1)+pm1(j,1)*0.5*cellsize(p)
                 pos(p,2)=pos(asubp(i),2)+pm1(j,2)*0.5*cellsize(p)
                 pos(p,3)=pos(asubp(i),3)+pm1(j,3)*0.5*cellsize(p)
                 bottom(p,1)=pos(p,1)-0.5*cellsize(p)
                 bottom(p,2)=pos(p,2)-0.5*cellsize(p)
                 bottom(p,3)=pos(p,3)-0.5*cellsize(p)
 100          CONTINUE

              nclist2=nclist2+nsubset

 110       CONTINUE

           nclist=nclist2
        
           DO 120 i=1,nclist
              celllist(i)=templist(i)
 120       CONTINUE

C   Find all subcells with one body; add bodies to tree.
C   ----------------------------------------------------
           DO 130 i=1,nbodlist
              templist(i)=ncells*(subindex(i)-1)+(parent(i)-nbodsmax)
              asubp(i)=subpvect(templist(i))
 130       CONTINUE

           CALL WHENEQ(nbodlist,asubp,1,1,isubset,nsubbod1)

CVD$ NODEPCHK
CDIR$ IVDEP
           DO 140 i=1,nsubbod1
              subpvect(templist(isubset(i)))=bodlist(isubset(i))
 140       CONTINUE

C   Place bodies in cells with more than one body on active list.
C   ------------------------------------------------------------

           CALL WHENIGT(nbodlist,asubp,1,1,isubset,nbodtemp)

           nbodlist=nbodtemp

           DO 150 i=1,nbodlist
              parent(i)=asubp(isubset(i))
              templist(i)=bodlist(isubset(i))
 150       CONTINUE

           DO 160 i=1,nbodlist
              bodlist(i)=templist(i)
 160       CONTINUE

           GO TO 200
        
        ENDIF 

        DO 230 i=1,nbodies
           bottom(i,1)=pos(i,1)
           bottom(i,2)=pos(i,2)
           bottom(i,3)=pos(i,3)
 230    CONTINUE

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE maketeps
C
C
C***********************************************************************
C
C
C     Subroutine to initialize the tree structure for computing
C     variable gravitational softening lengths for collisionless
C     particles.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

        LOGICAL firstc

        SAVE firstc

        DATA firstc/.TRUE./

C=======================================================================

        IF(firstc) THEN

           firstc=.FALSE.

           CALL setbox
C               ------
           CALL loadtree
C               --------
           CALL cellnumb
C               --------
        ENDIF

C   Compute coordinates of edges of cells.
C   --------------------------------------
        CALL celledge
C            --------

C   Determine which cells contain ingroup particles.
C   ------------------------------------------------
        CALL groupcel
C            --------

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE maketree
C
C
C***********************************************************************
C
C
C     Main routine to control initialization of the tree structure 
C     for computing the gravitational interaction.
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C=======================================================================

C   Set box properties.
C   -------------------
        CALL setbox
C            ------
 
C   Load bodies into the tree.
C   --------------------------
        CALL loadtree
C            --------

C   Compute properties of cells.
C   ----------------------------
        CALL hackcell
C            --------
 
        incellsg=incells
        ttree=tpos

        RETURN
        END
C***********************************************************************
C
C
                SUBROUTINE nearpeps(hsearch,pbody,npnear)
C
C
C***********************************************************************
C
C
C     Subroutine to search for nearest neighbors of body pbody within
C     hsearch softening lengths.  Vectorization is achieved by 
C     processing all cells at a given level in the tree simultaneously.
C     The neighbor list is passed back to the calling routine in the
C     vector nearlist, which is equivalenced to the common array
C     parent.
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER pbody,i,nsubdiv,nearlist(nbodsmax),npnear,nnodes,nkeep,
     &          k,ibody,jsubset(nbodsmax),njsubset,nodelist(ncells),
     &          keepnear(nbodsmax),keepstak(nbodsmax)
        LOGICAL testcrit
        REAL pbottom(ndim),ptop(ndim),hsearch,sradius,dx,dy,dz,
     &       xnode,ynode,znode

        EQUIVALENCE (nearlist(1),parent(1)),(jsubset(1),subindex(1)),
     &              (nodelist(1),celllist(1)),(keepnear(1),subindex(1)),
     &              (keepstak(1),asubp(1))

C=======================================================================

        npnear=0
        sradius=hsearch*epsvect(pbody)

        DO 10 k=1,3
           ptop(k)=pos(pbody,k)+sradius
           pbottom(k)=pos(pbody,k)-sradius
 10     CONTINUE

        nnodes=intone
        nodelist(1)=root

 20     CONTINUE

        IF(nnodes.GT.0) THEN

CVD$ NODEPCHK
CDIR$ IVDEP
           DO 30 i=1,nnodes
              xnode=bottom(nodelist(i),1)+0.5*cellsize(nodelist(i))
              ynode=bottom(nodelist(i),2)+0.5*cellsize(nodelist(i))
              znode=bottom(nodelist(i),3)+0.5*cellsize(nodelist(i))
              dx=0.0
              dy=0.0
              dz=0.0
              testcrit=     pbottom(1).LE.(bottom(nodelist(i),1)+
     &                                   cellsize(nodelist(i))+dx)
     &                 .AND.pbottom(2).LE.(bottom(nodelist(i),2)+
     &                                   cellsize(nodelist(i))+dy)
     &                 .AND.pbottom(3).LE.(bottom(nodelist(i),3)+
     &                                   cellsize(nodelist(i))+dz)
     &                 .AND.ptop(1).GE.(bottom(nodelist(i),1)+dx)
     &                 .AND.ptop(2).GE.(bottom(nodelist(i),2)+dy)
     &                 .AND.ptop(3).GE.(bottom(nodelist(i),3)+dz)
     &                 .AND.nodelist(i).NE.pbody
              IF(testcrit.AND.nodelist(i).LE.nbodsmax) THEN
                 keepnear(i)=2
              ELSE
                 keepnear(i)= -2
              ENDIF
              IF(testcrit.AND.nodelist(i).GT.nbodsmax) THEN
                 keepstak(i)=2
              ELSE
                 keepstak(i)= -2
              ENDIF
 30        CONTINUE

           CALL WHENIGT(nnodes,keepnear,1,0,isubset,nkeep)

           DO 35 i=1,nkeep
              ibody=nodelist(isubset(i))
              dx=pos(pbody,1)-pos(ibody,1)
              dy=pos(pbody,2)-pos(ibody,2)
              dz=pos(pbody,3)-pos(ibody,3)
              tempvect(i)=dx**2+dy**2+dz**2
 35        CONTINUE

           CALL WHENFLT(nkeep,tempvect,1,sradius**2,jsubset,njsubset)

           nkeep=njsubset

CVD$ NODEPCHK
CDIR$ IVDEP
           DO 40 i=1,nkeep
              nearlist(npnear+i)=nodelist(isubset(jsubset(i)))
 40        CONTINUE

           npnear=npnear+nkeep

           CALL WHENIGT(nnodes,keepstak,1,0,isubset,nsubdiv)

           IF(8*nsubdiv.GT.nbodsmax.OR.8*nsubdiv.GT.ncells)
     &        CALL terror(' asubp overflow in nearpeps ')
C                  ------

CVD$ NODEPCHK
CDIR$ IVDEP
           DO 50 i=1,nsubdiv
              asubp(i)=subp(nodelist(isubset(i)),1)
              asubp(i+nsubdiv)=subp(nodelist(isubset(i)),2)
              asubp(i+2*nsubdiv)=subp(nodelist(isubset(i)),3)
              asubp(i+3*nsubdiv)=subp(nodelist(isubset(i)),4)
              asubp(i+4*nsubdiv)=subp(nodelist(isubset(i)),5)
              asubp(i+5*nsubdiv)=subp(nodelist(isubset(i)),6)
              asubp(i+6*nsubdiv)=subp(nodelist(isubset(i)),7)
              asubp(i+7*nsubdiv)=subp(nodelist(isubset(i)),8)
 50        CONTINUE

           CALL WHENNE(8*nsubdiv,asubp,1,0,isubset,nnodes)

           DO 60 i=1,nnodes
              nodelist(i)=asubp(isubset(i))
 60        CONTINUE

           GO TO 20

        ENDIF

        RETURN
        END
C***********************************************************************
C
C
                         SUBROUTINE neighcol(pc)
C
C
C***********************************************************************
C
C
C     Subroutine to compute number of near neighbors for all
C     collisionless particles.  The list of neighbors is returned 
C     from the subroutine findnear in the vector nearlist, which is 
C     equivalenced to the common array bodlist.  If argument pc is
C     'predict' then this subroutine performs a simple neighbor
C     search for a given set of softening lengths.  If pc is 'correct'
C     this subroutine will additionally adjust softening lengths so
C     that the number of neighbors is approximately nsvolume.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        CHARACTER*7 pc
        INTEGER p,i,pbody,npnear,np,inear,jsubset(nbodsmax),
     &          nearlist(nbodsmax),keepnear(nbodsmax),
     &          nearpb(nbodsmax),nnearbod
        LOGICAL savecrit
        REAL epsinv,testnn,dx,dy,dz

        EQUIVALENCE (nearlist(1),bodlist(1)),(jsubset(1),subindex(1)),
     &              (keepnear(1),asubp(1)),(nearpb(1),parent(1))

C=======================================================================

C   Initialize neighbor diagnostics.
C   --------------------------------
        nstot=0
        nsmin=nbodies
        nsmax=0

C   Find nearest neighbors of grouped cells.
C   ----------------------------------------

        DO 100 p=1,ngroups

C   Find nearest neighbors.
C   -----------------------
           CALL findnear(p,npnear)
C               --------

           DO 70 np=pgroupb(p),pgroupb(p+1)-1

              pbody=groupbod(np)
              epsinv=1./epsvect(pbody)

              DO 20 i=1,npnear
                 dx=pos(pbody,1)-pos(nearlist(i),1)
                 dy=pos(pbody,2)-pos(nearlist(i),2)
                 dz=pos(pbody,3)-pos(nearlist(i),3)
                 tempvect(i)=(dx**2+dy**2+dz**2)*epsinv*epsinv
                 IF(nearlist(i).EQ.pbody) tempvect(i)=5.0
 20           CONTINUE

C   Filter out neighbors further than two smoothing lengths.
C   --------------------------------------------------------
              CALL WHENFLT(npnear,tempvect,1,four,isubset,inear)

              testnn=ABS(REAL(inear-nsvolume)/REAL(nsvolume))
              savecrit=(.NOT.variabls).OR.(testnn.LE.nsvtol).OR.
     &                 (pc.EQ.'predict')

              IF(savecrit) THEN

                 nstot=nstot+inear
                 IF(inear.LT.nsmin) nsmin=inear
                 IF(inear.GT.nsmax) nsmax=inear
                 nnear(pbody)=inear

              ELSE

                 IF(inear.GT.nsvolume) THEN

                    DO 50 i=1,inear
                       nearpb(i)=nearlist(isubset(i))
 50                 CONTINUE

                    CALL reducee(pbody,inear,nnearbod)
C                        -------

                 ELSE

                    CALL enlargee(pbody,nnearbod)
C                        --------
                 ENDIF

              ENDIF

 70        CONTINUE

 100    CONTINUE

        nsavg=nstot/nbodies

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE outbods
C
C
C***********************************************************************
C
C
C     Subroutine to output the body data.
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        CHARACTER*3 sstring
        CHARACTER*7 filename
        CHARACTER*8 filepar
        CHARACTER*10 nstring
        INTEGER p,ndimo,istring,nsnap,k

        SAVE nsnap,nstring

        DATA nsnap/0/,nstring/'0123456789'/

C=======================================================================
 
        nsnap=nsnap+1

        sstring(1:1)=nstring(1+nsnap/100:1+nsnap/100)
        istring=1+MOD(nsnap,100)/10
        sstring(2:2)=nstring(istring:istring)
        istring=1+MOD(nsnap,10)
        sstring(3:3)=nstring(istring:istring)
        filepar=obodfile
        filename=filepar(1:4)//sstring(1:3)

        OPEN(UNIT=ubodsout,FILE=filename,STATUS='NEW')

        ndimo=ndim

        WRITE(ubodsout,200) nbodies
        WRITE(ubodsout,200) ndimo
        WRITE(ubodsout,210) tnow
 
        DO 10 p=1,nbodies
           WRITE(ubodsout,210) mass(p)
 10     CONTINUE

        DO 20 p=1,nbodies
           WRITE(ubodsout,210) pos(p,1),pos(p,2),pos(p,3)
 20     CONTINUE


        DO 30 p=1,nbodies
           WRITE(ubodsout,210) vel(p,1),vel(p,2),vel(p,3)
 30     CONTINUE

        DO 40 p=1,nbodies
           WRITE(ubodsout,210) epsvect(p)
 40     CONTINUE

 200    FORMAT(1x,5(1i6))
 210    FORMAT(1x,10(1pe14.6))

        CLOSE(UNIT=ubodsout)

        RETURN
        END
C***********************************************************************
C
C
                         SUBROUTINE outcpu
C
C
C***********************************************************************
C
C
C     Subroutine to output cpu timing data to the log file.
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C=======================================================================
 
        OPEN(UNIT=ulog,FILE=logfile,STATUS='OLD',ACCESS='DIRECT',
     &       FORM='FORMATTED',RECL=81)

C   Output timing data to the log file.
C   -----------------------------------
        cputime=cputime1-cputime0

        ireclog=ireclog+1
        WRITE(ulog,10,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,10,REC=ireclog)

        ireclog=ireclog+1
        WRITE(ulog,20,REC=ireclog) cputime
 
        CLOSE(UNIT=ulog)

 10     FORMAT(' ')
 20     FORMAT(' Total cpu time used (seconds) : ',1pe12.4)
 
        RETURN
        END
C***********************************************************************
C
C
                         SUBROUTINE outdump(iopt)
C
C
C***********************************************************************
C
C
C     Subroutine to dump state of system to an ascii data file.  The
C     argument iopt indicates whether the data is output to a normal
C     SYSDUMP file (iopt=0) or a crash dump file (iopt=1).
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p,iopt,udumpout

C=======================================================================
 
        IF(iopt.EQ.0) THEN
           OPEN(UNIT=uboddump,FILE=dumpfile,STATUS='UNKNOWN')
           udumpout=uboddump
        ELSE
           OPEN(UNIT=ucrash,FILE=crashfil,STATUS='NEW')
           udumpout=ucrash
        ENDIF

C   Output system state.
C   --------------------

        WRITE(udumpout,500) nbodies

        WRITE(udumpout,510) tnow,tpos
        WRITE(udumpout,*) dtime
 
        WRITE(udumpout,*) mtot
        WRITE(udumpout,*) ektot
        WRITE(udumpout,*) eptot
        WRITE(udumpout,*) esofttot

        WRITE(udumpout,*) upbin
        WRITE(udumpout,*) stime
        WRITE(udumpout,*) tsteppos
        WRITE(udumpout,*) endstep

        DO 70 p=1,nbodies
           WRITE(udumpout,510) mass(p)
 70     CONTINUE

        DO 80 p=1,nbodies
           WRITE(udumpout,510) pos(p,1),pos(p,2),pos(p,3)
 80     CONTINUE

        DO 90 p=1,nbodies
           WRITE(udumpout,510) vel(p,1),vel(p,2),vel(p,3)
 90     CONTINUE

        DO 105 p=1,nbodies
           WRITE(udumpout,510) epsvect(p)
 105    CONTINUE

        DO 210 p=1,nbodies
           WRITE(udumpout,500) nnear(p)
 210    CONTINUE

        DO 220 p=1,nbodies
           WRITE(udumpout,510) phi(p)
 220    CONTINUE

        DO 230 p=1,nbodies
           WRITE(udumpout,510) acc(p,1),acc(p,2),acc(p,3)
 230    CONTINUE

        DO 240 p=1,nbodies
           WRITE(udumpout,500) itimestp(p)
 240    CONTINUE

        DO 250 p=1,nbodies
           WRITE(udumpout,500) otimestp(p)
 250    CONTINUE

        CLOSE(UNIT=udumpout)

 500    FORMAT(1x,10(1i10))
 510    FORMAT(1x,10(1pe22.14))

        RETURN
        END
C***********************************************************************
C
C
                  SUBROUTINE outenrgy(am,cmpos,cmvel)
C
C
C***********************************************************************
C
C
C     Subroutine to output diagnostic data to the log file.
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER k
        REAL am(3),cmpos(ndim),cmvel(ndim),cpunew,cpuold,cpustep,
     &       second

        SAVE cpuold

        DATA cpuold/0.0/

C=======================================================================

        OPEN(UNIT=ulog,FILE=logfile,STATUS='OLD',ACCESS='DIRECT',
     &       FORM='FORMATTED',RECL=81)

C-----------------------------------------------------------------------
C   Write mass, energy, angular momentum, center of mass quantities.
C-----------------------------------------------------------------------
 
        ireclog=ireclog+1
        WRITE(ulog,5,REC=ireclog)

        ireclog=ireclog+1
        WRITE(ulog,7,REC=ireclog) mtot

        ireclog=ireclog+1
        WRITE(ulog,20,REC=ireclog) etot,ektot,eptot
        ireclog=ireclog+1
        WRITE(ulog,22,REC=ireclog) esofttot

        ireclog=ireclog+1
        WRITE(ulog,40,REC=ireclog) am(1),am(2),am(3)
        ireclog=ireclog+1
        WRITE(ulog,50,REC=ireclog) (cmpos(k),k=1,ndim)
        ireclog=ireclog+1
        WRITE(ulog,60,REC=ireclog) (cmvel(k),k=1,ndim)

        cpunew=SECOND()
        cpustep=cpunew-cpuold
        cpuold=cpunew

        ireclog=ireclog+1
        WRITE(ulog,5,REC=ireclog)

        ireclog=ireclog+1
        WRITE(ulog,70,REC=ireclog) cpustep
 
        ireclog=ireclog+1
        WRITE(ulog,5,REC=ireclog)

        CLOSE(UNIT=ulog)

 5      FORMAT(' ')
 7      FORMAT(7x,'mtot = ',2(1pe12.4))
 20     FORMAT(7x,'e, ek, ep = ',3(1pe17.9))
 22     FORMAT(7x,'es = ',2(1pe17.9))
 40     FORMAT(7x,'amx, amy, amz = ',3(1pe17.9))
 50     FORMAT(7x,'cmpos = ',3(1pe17.9))
 60     FORMAT(7x,'cmvel = ',3(1pe17.9))
 70     FORMAT(10x,'cpu time per step = ',1pe12.4)
 
        RETURN
        END
C***********************************************************************
C
C
                       SUBROUTINE outerror(message)
C
C
C***********************************************************************
C
C
C     Subroutine to output error messages to the log file.
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        CHARACTER*(*) message

C=======================================================================

        OPEN(UNIT=ulog,FILE=logfile,STATUS='OLD',ACCESS='DIRECT',
     &       FORM='FORMATTED',RECL=81)

C   Write the message.
C   ------------------
        ireclog=ireclog+1
        WRITE(ulog,10,REC=ireclog)

        ireclog=ireclog+1
        WRITE(ulog,40,REC=ireclog)

        ireclog=ireclog+1
        WRITE(ulog,10,REC=ireclog)

        ireclog=ireclog+1
        WRITE(ulog,50,REC=ireclog) message

        ireclog=ireclog+1
        WRITE(ulog,10,REC=ireclog)

        ireclog=ireclog+1
        WRITE(ulog,40,REC=ireclog)
 
        CLOSE(UNIT=ulog)

 10     FORMAT(' ')
 40     FORMAT(1x,72('*'))
 50     FORMAT(a)

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE outhead
C
C
C***********************************************************************
C
C
C     Subroutine to output a standard header to the log file.
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C=======================================================================
 
        OPEN(UNIT=ulog,FILE=logfile,STATUS='OLD',ACCESS='DIRECT',
     &       FORM='FORMATTED',RECL=81)

        ireclog=ireclog+1
        WRITE(ulog,10,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,20,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,20,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,25,REC=ireclog) headline
        ireclog=ireclog+1
        WRITE(ulog,20,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,10,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,20,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,28,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,29,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,30,REC=ireclog) nbodies,nsteps,noutbod,noutlog
        ireclog=ireclog+1
        WRITE(ulog,20,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,40,REC=ireclog) dtime,eps,usequad,tol
        ireclog=ireclog+1
        WRITE(ulog,20,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,42,REC=ireclog) dgravsum,variabls,nsvolume
        ireclog=ireclog+1
        WRITE(ulog,20,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,70,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,80,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,20,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,90,REC=ireclog) inittbin,mintstep,etol,ntvector
        ireclog=ireclog+1
        WRITE(ulog,20,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,10,REC=ireclog)
 
        CLOSE(UNIT=ulog)

 10     FORMAT(1x,72('*'))
 20     FORMAT(1x,'*',70(' '),'*')
 25     FORMAT(1x,'*',10x,1a50,10x,'*')
 28     FORMAT(1x,'*',4x,'Input parameters:',49x,'*')
 29     FORMAT(1x,'*',4x,'----------------',50x,'*')
 30     FORMAT(1x,'*',8x,'nbodies=',1i7,2x,'nsteps=',1i4,4x,
     &         'noutbods=',1i4,2x,'noutlog=',1i4,3x,'*')
 40     FORMAT(1x,'*',8x,'dtime=',1pe10.3,1x,'eps=',1pe10.3,1x,
     &         'usequad=',1l1,6x,'tol=',0pf5.2,6x,'*')
 42     FORMAT(1x,'*',8x,'dgravsum=',1l1,2x,'variabls=',1l1,2x,
     &         'nsvolume=',1i7,22x,'*')
 70     FORMAT(1x,'*',4x,'Time step parameters:',45x,'*')
 80     FORMAT(1x,'*',4x,'--------------------',46x,'*')
 90     FORMAT(1x,'*',8x,'inittbin=',1i5,2x,'mintstep=',1i5,2x,
     &         'etol=',1f7.3,2x,'ntvector=',1i5,2x,'*')

        RETURN
        END
C***********************************************************************
C
C
                        SUBROUTINE outlog(istep)
C
C
C***********************************************************************
C
C
C     Subroutine to monitor status of the program by writing to
C     the log file.
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER istep,i,tstepbin(20),p
        REAL tistep

C=======================================================================
 
C   If first call, write header.
C   ----------------------------
        IF(istep.EQ.0) CALL outhead
C                           -------
 
C-----------------------------------------------------------------------
C   Output system time and force evaluation diagnostics.
C-----------------------------------------------------------------------

        OPEN(UNIT=ulog,FILE=logfile,STATUS='OLD',ACCESS='DIRECT',
     &       FORM='FORMATTED',RECL=81)
 
        ireclog=ireclog+1
        WRITE(ulog,65,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,65,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,75,REC=ireclog) tnow,incellsg
        ireclog=ireclog+1
        WRITE(ulog,65,REC=ireclog)

        ireclog=ireclog+1
        WRITE(ulog,80,REC=ireclog) nttot,ntmin,ntmax,ntavg
        ireclog=ireclog+1
        WRITE(ulog,90,REC=ireclog) nntot,nnmin,nnmax,nnavg

        IF(variabls) THEN
           ireclog=ireclog+1
           WRITE(ulog,95,REC=ireclog) nstot,nsmin,nsmax,nsavg
        ENDIF

        ireclog=ireclog+1 
        WRITE(ulog,65,REC=ireclog)

        DO 40 i=1,20
           tstepbin(i)=0
 40     CONTINUE

        DO 50 p=1,nbodies
           tistep=itimestp(p)
           tistep=0.5+LOG(tistep)/log2
           templist(p)=INT(tistep)+1
 50     CONTINUE

        DO 60 p=1,nbodies
           tstepbin(templist(p))=tstepbin(templist(p))+1
 60     CONTINUE

        ireclog=ireclog+1
        WRITE(ulog,65,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,120,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,65,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,110,REC=ireclog) (tstepbin(i),i=1,8)
        ireclog=ireclog+1
        WRITE(ulog,110,REC=ireclog) (tstepbin(i),i=9,16)

        CLOSE(UNIT=ulog)

 65     FORMAT(' ')
 70     FORMAT(2x,'time: ',1pe11.3,5x,'aexp: ',1pe11.3,5x,'Z: ',
     &         0pf5.2,5x,'ncells: ',1i5)
 75     FORMAT(2x,'time: ',1pe12.4,5x,'ncells: ',1i5)
 80     FORMAT(7x,'nttot, min, max, avg = ',1i8,5x,1i5,5x,1i5,5x,1i5)
 90     FORMAT(7x,'nntot, min, max, avg = ',1i8,5x,1i5,5x,1i5,5x,1i5)
 95     FORMAT(7x,'nstot, min, max, avg = ',1i8,5x,1i5,5x,1i5,5x,1i5)
 110    FORMAT(12x,8i6)
 120    FORMAT(7x,'distribution of time steps : ')

        RETURN
        END
C***********************************************************************
C
C
                        SUBROUTINE outstate(n)
C
C
C***********************************************************************
C
C
C     Subroutine to output information about the system state to
C     the log and body data files.
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER n,ioutcosm

        SAVE ioutcosm

        DATA ioutcosm/1/

C=======================================================================

        CALL outterm(' step completed: ',n)
C            -------

        IF(n.EQ.0) THEN

           CALL outlog(0)
C               ------
           CALL energy
C               ------
           CALL outbods
C               -------

        ELSE


c-------------------------
c         animation

           IF(MOD(n,noutani).EQ.0) 
     +          CALL anitul(n,tnow,nbodies,pos)
           

           IF(MOD(n,noutlog).EQ.0.OR.(MOD(n,noutbod).EQ.0)) THEN

              IF(MOD(n,noutbod).EQ.0) CALL outdump(0)
C                                          -------
              CALL corrpos('correct')
C                  -------
              CALL zeropot
C                  -------
              CALL gravity('pot ')
C                  -------

              IF(MOD(n,noutlog).EQ.0) THEN

                 CALL outlog(n)
C                     ------
                 CALL energy
C                     ------
              ENDIF

              IF(MOD(n,noutbod).EQ.0) THEN 
                 CALL outbods
C                     -------
              ENDIF

              CALL corrpos('reset  ')
C                  -------
           ENDIF

        ENDIF
 
        RETURN
        END
C***********************************************************************
C
C
                      SUBROUTINE outterm(message,n)
C
C
C***********************************************************************
C
C
C     Subroutine to output a message to the terminal and to the
C     terminal emulation file.
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        CHARACTER*(*) message
        INTEGER n

C=======================================================================
 
        OPEN(UNIT=utermfil,FILE=termfile,STATUS='OLD')

C   Write the message.
C   ------------------

        IF(n.GE.0) THEN
           WRITE(uterm,*) message,n
           WRITE(utermfil,*) message,n
        ELSE
           WRITE(uterm,40)
           WRITE(uterm,50) message 
           WRITE(uterm,40)
           WRITE(utermfil,40)
           WRITE(utermfil,50) message 
           WRITE(utermfil,40)
        ENDIF

 40     FORMAT(/,1x,72('*'))
 50     FORMAT(/,a)

        CLOSE(UNIT=utermfil)

        RETURN
        END
C***********************************************************************
C
C
               SUBROUTINE reducee(pbody,inear,nnearbod)
C
C
C***********************************************************************
C
C
C     Subroutine to decrease softening length of body pbody so that
C     the softening volume contains roughly nsvolume neighbors.  The
C     list of nearest neighbors found in stepeps is passed through
C     the vector nearpb, which is equivalenced to the common array
C     parent.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER i,pbody,inear,nhtest,nearpb(nbodsmax),nnearbod,
     &          jsubset(nbodsmax),keepnear(nbodsmax)
        REAL hmin,hmax,htest,fourhsm2,dx,dy,dz

        EQUIVALENCE (nearpb(1),parent(1)),(jsubset(1),subindex(1)),
     &              (keepnear(1),asubp(1))

C=======================================================================

        hmax=epsvect(pbody)

        DO 5 i=1,inear
           isubset(i)=i
           dx=pos(pbody,1)-pos(nearpb(i),1)
           dy=pos(pbody,2)-pos(nearpb(i),2)
           dz=pos(pbody,3)-pos(nearpb(i),3)
           tempvect(i)=dx**2+dy**2+dz**2
 5      CONTINUE

        hmin=0.5*hmax

 7      CONTINUE

        htest=hmin

        fourhsm2=four*htest*htest

        CALL WHENFLT(inear,tempvect,1,fourhsm2,isubset,nhtest)

        IF(nhtest.GE.nsvolume) THEN
           hmin=0.5*hmin
           GO TO 7
        ENDIF

 10     CONTINUE

        IF(ABS(REAL(nhtest-nsvolume)).GT.nsvtol*REAL(nsvolume)) THEN

           htest=0.5*(hmin+hmax)
           fourhsm2=four*htest*htest

           CALL WHENFLT(inear,tempvect,1,fourhsm2,isubset,nhtest)

           IF(nhtest.GT.nsvolume) THEN
              hmax=htest
           ELSE
              hmin=htest
           ENDIF

           GO TO 10

        ENDIF

        nstot=nstot+nhtest
        IF(nhtest.LT.nsmin) nsmin=nhtest
        IF(nhtest.GT.nsmax) nsmax=nhtest
        nnear(pbody)=nhtest
        epsvect(pbody)=htest

        RETURN
        END
C***********************************************************************
C
C
                           SUBROUTINE setbox
C
C
C***********************************************************************
C
C
C     Subroutine to adjust system box so that it contains all bodies.
C     The local variable rebox indicates whether a resizing of the 
C     system box is to take place.  The variables posmin and posmax 
C     are the minimum and maximum coordinates of bodies in each 
C     dimension.  
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER k,ismin,ismax
        LOGICAL rebox
        REAL posmin(ndim),posmax(ndim),posx(nbodsmax),posy(nbodsmax),
     &       posz(nbodsmax)

        EQUIVALENCE (posx(1),pos(1,1)),(posy(1),pos(1,2)),
     &              (posz(1),pos(1,3))

        SAVE rebox

        DATA rebox/.TRUE./

C=======================================================================

C   Determine minimum and maximum coordinates of bodies.
C   ----------------------------------------------------

        posmin(1)=posx(ISMIN(nbodies,posx,1))
        posmin(2)=posy(ISMIN(nbodies,posy,1))
        posmin(3)=posz(ISMIN(nbodies,posz,1))
        posmax(1)=posx(ISMAX(nbodies,posx,1))
        posmax(2)=posy(ISMAX(nbodies,posy,1))
        posmax(3)=posz(ISMAX(nbodies,posz,1))

C   Determine if a resizing is required.
C   ------------------------------------
        DO 50 k=1,ndim
           IF(rmin(k).GT.posmin(k).OR.rmin(k)+rsize.LT.posmax(k)) 
     &        rebox=.TRUE.
 50     CONTINUE

C   If a resizing is necessary, recompute rsize and rmin.
C   -----------------------------------------------------

        IF(rebox) THEN

           DO 70 k=1,ndim
              IF(rsize.LT.posmax(k)-posmin(k)) rsize=posmax(k)-posmin(k)
 70        CONTINUE

           DO 80 k=1,ndim
              rmin(k)=0.5*(posmin(k)+posmax(k))-0.5*rsize
 80        CONTINUE

        ENDIF

        rebox=.FALSE.

        RETURN
        END
C***********************************************************************
C
C
                SUBROUTINE srchbox(p,npnear,pbottom,ptop)
C
C
C***********************************************************************
C
C
C     Subroutine to initialize the search box for the nearest
C     neighbor detection for a group of particles specified by
C     the argument p.  The list of neighbors is returned through
C     the vector nearlist, which is equivalenced to the common
C     array bodlist.
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p,j,k,nearlist(nbodsmax),npgroup,ismin,ismax,npnear
        REAL pbottom(ndim),ptop(ndim),rsearch

        EQUIVALENCE (nearlist(1),bodlist(1))

C=======================================================================
 
        npgroup=pgroupb(p+1)-pgroupb(p)

        DO 30 k=1,3

CVD$ NODEPCHK
CDIR$ IVDEP
           DO 10 j=pgroupb(p),pgroupb(p+1)-1
              rsearch=epsvect(groupbod(j))
              tempvect(j-pgroupb(p)+1)=pos(groupbod(j),k)-2.*rsearch
 10        CONTINUE

           pbottom(k)=tempvect(ISMIN(npgroup,tempvect,1))

CVD$ NODEPCHK
CDIR$ IVDEP
           DO 20 j=pgroupb(p),pgroupb(p+1)-1
              rsearch=epsvect(groupbod(j))
              tempvect(j-pgroupb(p)+1)=pos(groupbod(j),k)+2.*rsearch
 20        CONTINUE

           ptop(k)=tempvect(ISMAX(npgroup,tempvect,1))

 30     CONTINUE

CVD$ NODEPCHK
CDIR$ IVDEP
        DO 40 j=pgroupb(p),pgroupb(p+1)-1
           nearlist(npnear+j-pgroupb(p)+1)=groupbod(j)
 40     CONTINUE

        npnear=npnear+npgroup

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE startout
C
C
C***********************************************************************
C
C
C     Subroutine to open disk files for subsequent input/output.
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C=======================================================================
 
C   Open log file.
C   --------------
        OPEN(UNIT=ulog,FILE=logfile,STATUS='NEW',ACCESS='DIRECT',
     &       FORM='FORMATTED',RECL=81)

        ireclog=0

        WRITE(ulog,10,REC=1)
 10     FORMAT(' Start of logfile output ')

        CLOSE(UNIT=ulog)
 
C   Create terminal emulation file.
C   -------------------------------
        OPEN(UNIT=utermfil,FILE=termfile,STATUS='UNKNOWN')
        WRITE(utermfil,*) ' Start of output '
        CLOSE(UNIT=utermfil)

        RETURN
        END
C***********************************************************************
C
C
                           SUBROUTINE stepeps
C
C
C***********************************************************************
C
C
C     Subroutine to update softening lengths for all collisionless
C     particles.  The list of neighbors is returned from the 
C     subroutine findnear in the vector nearlist, which is 
C     equivalenced to the common array bodlist.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p,i,pbody,npnear,inear,nnearbod,nearpb(nbodsmax),
     &          np,jsubset(nbodsmax),nearlist(nbodsmax),
     &          keepnear(nbodsmax),nnearfix
        REAL testnn,foureps2,dx,dy,dz

        EQUIVALENCE (nearlist(1),bodlist(1)),(jsubset(1),subindex(1)),
     &              (nearpb(1),parent(1)),(keepnear(1),asubp(1))

C=======================================================================

C   Extrapolate softening lengths.
C   ------------------------------
        DO 10 p=1,nbodies
           IF(nnear(p).EQ.0) THEN
              nnearfix=1
           ELSE
              nnearfix=0
           ENDIF
           epsvect(p)=epsvect(p)*0.5*((REAL(nsvolume)/REAL(nnear(p)+
     &                nnearfix))**(one/3.)+one)
 10     CONTINUE

C   Initialize neighbor diagnostics.
C   --------------------------------
        nstot=0
        nsmin=nbodies
        nsmax=0

C   Find nearest neighbors of grouped cells.
C   ----------------------------------------

        DO 100 p=1,ngroups

C   Find nearest neighbors.
C   -----------------------
           CALL findnear(p,npnear)
C               --------

           DO 70 np=pgroupb(p),pgroupb(p+1)-1

              pbody=groupbod(np)

              foureps2=four*epsvect(pbody)*epsvect(pbody)

              DO 20 i=1,npnear
                 dx=pos(pbody,1)-pos(nearlist(i),1)
                 dy=pos(pbody,2)-pos(nearlist(i),2)
                 dz=pos(pbody,3)-pos(nearlist(i),3)
                 tempvect(i)=dx**2+dy**2+dz**2
                 IF(nearlist(i).EQ.pbody) tempvect(i)=foureps2+one
 20           CONTINUE

C   Filter out neighbors further than two smoothing lengths.
C   --------------------------------------------------------

              CALL WHENFLT(npnear,tempvect,1,foureps2,isubset,inear)

              testnn=ABS(REAL(inear-nsvolume)/REAL(nsvolume))

              IF(testnn.LE.nsvtol) THEN

                 nstot=nstot+inear
                 IF(inear.LT.nsmin) nsmin=inear
                 IF(inear.GT.nsmax) nsmax=inear
                 nnear(pbody)=inear

              ELSE

                 IF(inear.GT.nsvolume) THEN

                    DO 50 i=1,inear
                       nearpb(i)=nearlist(isubset(i))
 50                 CONTINUE
 
                    CALL reducee(pbody,inear,nnearbod)
C                        -------

                 ELSE

                    CALL enlargee(pbody,nnearbod)
C                        --------
                 ENDIF

              ENDIF

 70        CONTINUE

 100    CONTINUE

        nsavg=nstot/nbodies

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE steppos
C
C
C***********************************************************************
C
C
C     Subroutine to advance the positions of the bodies for a
C     timestep dt.  The local variable p is a pointer to the bodies.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p,k

C=======================================================================

C   Loop over all spatial coordinates for all bodies.
C   -------------------------------------------------

        DO 200 k=1,ndim
           DO 100 p=1,nbodies
              pos(p,k)=pos(p,k)+vel(p,k)*tsteppos
 100       CONTINUE
 200    CONTINUE

C   Update position time, system time.
C   ----------------------------------
        tpos=tpos+tsteppos
        tnow=tpos

        RETURN
        END
C***********************************************************************
C
C
                         SUBROUTINE stepsys(n)
C
C
C***********************************************************************
C
C
C     Subroutine to advance the state of the system by one large
C     timestep.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER n

C=======================================================================

 50     CONTINUE

           CALL timestep
C               --------

C   Update positions.
C   -----------------
           CALL steppos
C               -------

           IF(npactive.GT.0) THEN

C   Zero out acceleration, compute gravitational acceleration.
C   ----------------------------------------------------------
              IF(.NOT.endstep) THEN

                 CALL zeroacc
C                     -------
                 CALL gravity('acc ')
C                     -------
              ENDIF

           ENDIF
 
           IF(.NOT.endstep) THEN

              CALL stepvel
C                  -------
              GO TO 50

           ENDIF

C   Output system state.
C   --------------------
        CALL outstate(n)
C            --------

        RETURN
        END
C***********************************************************************
C
C
                         SUBROUTINE stepvel
C
C
C***********************************************************************
C
C
C     Subroutine to advance the velocities of the bodies for a
C     timestep dt.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p,k,i

C=======================================================================

C   Loop over all velocity components for all bodies.
C   -------------------------------------------------

        DO 200 k=1,ndim

CVD$ NODEPCHK
CDIR$ IVDEP
           DO 100 i=1,npactive
              p=pactive(i)
              vel(p,k)=vel(p,k)+acc(p,k)*dtime/itimestp(p)
 100       CONTINUE

 200    CONTINUE

        RETURN
        END
C***********************************************************************
C
C
                           SUBROUTINE stopout
C
C
C***********************************************************************
C
C
C     Subroutine to close the open output files.
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C=======================================================================
 
C   Close the open files.
C   ---------------------
C       CLOSE(UNIT=ulog)
 
        RETURN
        END
C***********************************************************************
C
C
                       SUBROUTINE terror(message)
C
C
C***********************************************************************
C
C
C     Subroutine to terminate the program as the result of a fatal
C     error, close the output files, and dump timing information.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        CHARACTER*(*) message
        INTEGER ierror
        REAL second

C=======================================================================

C   Write error message to the log file and to the terminal.
C   --------------------------------------------------------
        CALL outerror(message)
C            --------
        ierror=-1

        CALL outterm(message,ierror)
C            -------

        CALL outdump(1)
C            -------
        CALL outbods
C            -------

C-----------------------------------------------------------------------
C   Stop timing, output timing data, close files, terminate the
C   simulation.
C-----------------------------------------------------------------------

        cputime1=SECOND()

        CALL outcpu
C            ------
        CALL stopout
C            -------

        STOP
        END
C***********************************************************************
C
C
                          SUBROUTINE timestep
C
C
C***********************************************************************
C
C
C     Subroutine to find appropriate time step for each particle
C     for the variable time step scheme.  Also finds the smallest
C     timestep and determines which particles are to be moved during
C     the next time step
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER i,p,ttemp,ideal,nrem,npacts,npactl,tmin,upbint,ismax,
     &          nchange,idele
        LOGICAL testcrit
        REAL epart,vt,at,c1,tflag,epartp

C=======================================================================

        IF(ektot.LT.-eptot/6.) THEN
           epart= -eptot/(6.*mtot)
        ELSE
           epart=ektot/mtot
        ENDIF

        CALL WHENIGT(nbodies,itimestp,1,upbin,pactive,npactive)

        DO 7 i=1,npactive
           templist(i)=0
 7      CONTINUE

CVD$ NODEPCHK
CDIR$ IVDEP
        DO 10 i=1,npactive
           p=pactive(i)
           vt=vel(p,1)*vel(p,1)+vel(p,2)*vel(p,2)+vel(p,3)*vel(p,3)
           at=acc(p,1)*acc(p,1)+acc(p,2)*acc(p,2)+acc(p,3)*acc(p,3)
           epartp=0.5*vt+phi(p)
           epartp=ABS(epartp)
           epart=epartp
           idele=dtime*SQRT(at*vt)/(etol*epart)+1.
           IF(templist(i).LT.idele) templist(i)=idele
           testcrit=itimestp(p).LT.templist(i).AND.itimestp(p)
     &              .LT.mintstep
           IF(testcrit) THEN
              tflag=two
           ELSE
              tflag=zero
           ENDIF
           testcrit=(itimestp(p).NE.1).AND.(itimestp(p).GT.2*upbin)
     &              .AND.(itimestp(p)/2.GE.templist(i))
           IF(testcrit) THEN
              tempvect(i)=four
           ELSE
              tempvect(i)=zero
           ENDIF
           tempvect(i)=tempvect(i)+tflag
 10     CONTINUE

        CALL WHENFGT(npactive,tempvect,1,three,isubset,nchange)

CVD$ NODEPCHK
CDIR$ IVDEP
        DO 15 i=1,nchange
           p=pactive(isubset(i))
           itimestp(p)=itimestp(p)/2
           tempvect(isubset(i))=0.
 15     CONTINUE

        CALL WHENFGT(npactive,tempvect,1,one,isubset,nchange)

CVD$ NODEPCHK
CDIR$ IVDEP
        DO 16 i=1,nchange
           p=pactive(isubset(i))
           c1=templist(isubset(i))
           ideal=LOG(c1)/log2 + 1.
           ideal=2**ideal + .5
           IF(mintstep.LT.ideal) THEN
              itimestp(p)=mintstep
           ELSE
              itimestp(p)=ideal
           ENDIF
 16     CONTINUE

        DO 17 p=1,nbodies
           tempvect(p)=itimestp(p)
 17     CONTINUE

        tmin=itimestp(ISMAX(nbodies,tempvect,1))

        stime=stime + 1./(tmin*2.)
        upbin=tmin
        ttemp=2.*tmin*stime+.5

        IF(ttemp.GE.2*tmin) THEN

           tsteppos=dtime*(1.-stime+1./(tmin*2.))
           upbin=0
           stime=0.0
           npactive=nbodies
           endstep=.TRUE.

           DO 55 i=1,npactive
              pactive(i)=i
 55        CONTINUE

           RETURN
   
        ENDIF

 20     CONTINUE

        IF(upbin.GT.1.AND.MOD(ttemp,2).EQ.0) THEN
           upbin=upbin/2
           ttemp=ttemp/2
           GO TO 20
        ENDIF

        CALL WHENEQ(nbodies,itimestp,1,upbin,pactive,npactive)

        IF(npactive.GT.0) THEN

        upbint=upbin

 23        CONTINUE

           nrem=MOD(npactive,ntvector)
           nrem=ntvector-nrem
           npacts=MOD(nrem,ntvector)

           IF(npacts.NE.0.AND.upbint.NE.1) THEN
              npactl=npactive
              upbint=upbint/2
              IF(upbint.EQ.0) CALL terror('error in timestep')
C                                  ------

              CALL WHENEQ(nbodies,itimestp,1,upbint,bodlist,nchange)

              IF(nchange.LT.npacts) npacts=nchange
              npactive=npactive+npacts

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 30 i=1,npacts
                 p=bodlist(i)
                 pactive(i+npactl)=p
                 itimestp(p)=upbin
 30           CONTINUE
              GO TO 23
           ENDIF

CVD$ NODEPCHK
CDIR$ IVDEP
           DO 40 i=1,npactive
              p=pactive(i)
              isubset(i)=itimestp(p)/otimestp(p)
              otimestp(p)=itimestp(p)
 40        CONTINUE

           CALL WHENEQ(npactive,isubset,1,0,templist,nchange)

CVD$ NODEPCHK
CDIR$ IVDEP
           DO 45 i=1,nchange
              p=pactive(templist(i))
              c1=3./32.
              c1=c1*dtime*dtime/(itimestp(p)*itimestp(p))
              pos(p,1)=pos(p,1)+c1*acc(p,1)
              pos(p,2)=pos(p,2)+c1*acc(p,2)
              pos(p,3)=pos(p,3)+c1*acc(p,3)
 45        CONTINUE

           CALL WHENIGT(npactive,isubset,1,1,templist,nchange)

CVD$ NODEPCHK
CDIR$ IVDEP
           DO 50 i=1,nchange
              p=pactive(templist(i))
              c1=-1./8.*(1.-1./isubset(templist(i)))*
     &         (1.+1./isubset(templist(i)))*(isubset(templist(i))**2)
              c1=c1*dtime*dtime/(itimestp(p)*itimestp(p))
              pos(p,1)=pos(p,1)+c1*acc(p,1)
              pos(p,2)=pos(p,2)+c1*acc(p,2)
              pos(p,3)=pos(p,3)+c1*acc(p,3)
 50        CONTINUE

        ENDIF

        tsteppos=dtime/(2.*tmin)
        endstep=.FALSE.

        RETURN
        END
C***********************************************************************
C
C
                     SUBROUTINE treewalk(p,nterms)
C
C
C***********************************************************************
C
C
C     Subroutine to walk through the tree and accumulate the list of
C     interactions for body p.  The interaction list is passed back
C     to the calling subroutine accgrav in the vector iterms, which
C     is equivalenced to the common array bodlist.
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p,i,nnodes,nkeep,nsubdiv,nterms,iterms(nbodsmax),
     &          nodelist(ncells),keepterm(nbodsmax)
        LOGICAL tolcrit

        EQUIVALENCE (iterms(1),bodlist(1)),(nodelist(1),celllist(1)),
     &              (keepterm(1),parent(1))

C=======================================================================

C   Initialize list of cells to examine.
C   ------------------------------------
        nterms=0
        nnodes=intone
        nodelist(1)=root
     
 10     CONTINUE

C   Loop until no cells are left to examine.
C   ----------------------------------------
        IF(nnodes.GT.0) THEN

C   Apply tolerance criterion to list of cells.
C   -------------------------------------------

CVD$ NODEPCHK
CDIR$ IVDEP
           DO 20 i=1,nnodes
              tolcrit=(tol2*((pos(p,1)-pos(nodelist(i),1))**2+
     &                (pos(p,2)-pos(nodelist(i),2))**2+
     &                (pos(p,3)-pos(nodelist(i),3))**2)).GE.
     &                 cellsize(nodelist(i))**2
              IF(tolcrit) THEN
                 keepterm(i)=2
              ELSE
                 keepterm(i) = -2
              ENDIF
 20        CONTINUE

C-----------------------------------------------------------------------
C   Add cells which satisfy criterion to interaction list.  Note that,
C   depending on theta, self-interaction term may be included.
C-----------------------------------------------------------------------

           CALL WHENIGT(nnodes,keepterm,1,0,isubset,nkeep)

           IF(nterms+nkeep.GT.nbodsmax) 
     &        CALL terror(' array overflow in treewalk ')
C                  ------

CVD$ NODEPCHK
CDIR$ IVDEP
           DO 30 i=1,nkeep
              iterms(nterms+i)=nodelist(isubset(i))
 30        CONTINUE

           nterms=nterms+nkeep

C-----------------------------------------------------------------------
C   Add subcells of cells which fail tolerance criterion to list of
C   cells to examine.
C-----------------------------------------------------------------------

           CALL WHENILT(nnodes,keepterm,1,0,isubset,nsubdiv)

           IF(8*nsubdiv.GT.nbodsmax.OR.8*nsubdiv.GT.ncells)
     &        CALL terror(' asubp overflow in treewalk ')
C                  ------

CVD$ NODEPCHK
CDIR$ IVDEP
           DO 40 i=1,nsubdiv
              asubp(i)=subp(nodelist(isubset(i)),1)
              asubp(i+nsubdiv)=subp(nodelist(isubset(i)),2)
              asubp(i+2*nsubdiv)=subp(nodelist(isubset(i)),3)
              asubp(i+3*nsubdiv)=subp(nodelist(isubset(i)),4)
              asubp(i+4*nsubdiv)=subp(nodelist(isubset(i)),5)
              asubp(i+5*nsubdiv)=subp(nodelist(isubset(i)),6)
              asubp(i+6*nsubdiv)=subp(nodelist(isubset(i)),7)
              asubp(i+7*nsubdiv)=subp(nodelist(isubset(i)),8)
 40        CONTINUE

           IF(nsubdiv.GT.0) THEN
              CALL WHENNE(8*nsubdiv,asubp,1,0,isubset,nnodes)
           ELSE
              nnodes=intzero
           ENDIF

CVD$ NODEPCHK
CDIR$ IVDEP
           DO 60 i=1,nnodes
              nodelist(i)=asubp(isubset(i))
 60        CONTINUE

           GO TO 10

        ENDIF

        RETURN
        END
C***********************************************************************
C
C
                         SUBROUTINE zeroacc
C
C
C***********************************************************************
C
C
C     Subroutine to zero out acceleration.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p,i

C=======================================================================

CVD$ NODEPCHK
CDIR$ IVDEP
        DO 10 i=1,npactive
           p=pactive(i)
           acc(p,1)=0.
           acc(p,2)=0.
           acc(p,3)=0.
 10     CONTINUE

        RETURN
        END
C***********************************************************************
C
C
                         SUBROUTINE zeropot
C
C
C***********************************************************************
C
C
C     Subroutine to zero out potential.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p,i

C=======================================================================

CVD$ NODEPCHK
CDIR$ IVDEP
        DO 10 i=1,npactive
           p=pactive(i)
           phi(p)=0.
 10     CONTINUE

        IF(npactive.EQ.nbodies) esofttot=0.0

        RETURN
        END

