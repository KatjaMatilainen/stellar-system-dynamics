C=======================================================================
C
C
C                        INCLUDE FILE treedefs.h                       
C
C
C=======================================================================
C
C
C     Parameter declarations, allocation of array storage, common
C     block definitions.
C
C
C=======================================================================

        CHARACTER*50 headline
        INTEGER root,subp,nbodies,incells,nttot,ntmin,ntmax,ntavg,
     &          nsteps,noutbod,noutlog,ndim,nsubcell,nbodsmax,ncells,
     &          nbodcell,nbods1,incellsg,intzero,intone,inttwo,
     &          intthree,intfour,nsvolume,nstot,nsmin,nsmax,nsavg,
     &          ninterp
        LOGICAL usequad,selfgrav,variabls,dgravsum,inpteps,restart
        REAL mass,tol,tol2,eps,rsize,rmin,phi,pos,vel,acc,quad,tnow,
     &       tpos,dtime,dtime2,tiny,one,two,four,cputime0,cputime1,etot,
     &       cputime,zero,zero02,cellsize,minusone,ttree,mtot,ektot,
     &       eptot,three,log2,minustwo,epsvect,onehalf,nsvtol,esofttot,
     &       phsmooth,acsmooth

        PARAMETER(ndim=3,nsubcell=2**ndim)
        PARAMETER(nbodsmax=36864,ncells=20000,ninterp=30000)
        PARAMETER(nbodcell=nbodsmax+ncells,nbods1=nbodsmax+1)

        COMMON/paramcom/nbodies,tol,tol2,eps,usequad,selfgrav,variabls,
     &                  nsvolume,nsvtol,dgravsum,inpteps,restart
        COMMON/msgcom/headline
        COMMON/cellcom/rsize,rmin(ndim),incells,incellsg
        COMMON/pointers/root,subp(nbods1:nbodcell,1:nsubcell)
        COMMON/bodycell/mass(1:nbodcell),phi(1:nbodsmax),
     &                  pos(1:nbodcell,1:ndim),cellsize(1:nbodcell),
     &                  vel(1:nbodsmax,1:ndim),acc(1:nbodsmax,1:ndim),
     &                  epsvect(1:nbodcell)
        COMMON/quadcom/quad(nbods1:nbodcell,1:2*ndim-1)
        COMMON/forcecom/nttot,ntmin,ntmax,ntavg
        COMMON/softcom/nstot,nsmin,nsmax,nsavg
        COMMON/timecom/nsteps,noutbod,noutlog,tnow,tpos,dtime,dtime2,
     &                 ttree
        COMMON/cpucom/cputime0,cputime1,cputime
        COMMON/misccom/tiny,zero,one,two,four,minusone,zero02,three,
     &                 log2,minustwo,intone,inttwo,intthree,intfour,
     &                 intzero,onehalf
        COMMON/enrgycom/mtot,etot,ektot,eptot,esofttot
        COMMON/gravcom/phsmooth(0:1+ninterp),acsmooth(0:1+ninterp)

C=======================================================================
C   Definitions specific to input/output.
C=======================================================================
        INTEGER uterm,upars,ulog,ubodsin,ubodsout,uboddump,utermfil,
     &          ireclog,ucrash,uindump
        CHARACTER*8 parsfile,logfile,ibodfile,obodfile,dumpfile,
     &              termfile,crashfil,indumpf

        PARAMETER(uterm=6,upars=10,ulog=11,ubodsin=12,ubodsout=13,
     &            uboddump=14,utermfil=15,ucrash=16,uindump=17)
        PARAMETER(parsfile='TREEPAR',logfile='TREELOG',
     &            ibodfile='TREEBI',obodfile='SNAPxxxx',
     &            dumpfile='SYSDUMP',termfile='TREEOUT',
     &            crashfil='CDUMP',indumpf='DUMPIN')

        COMMON/iocommon/ireclog

C=======================================================================
C   Definitions specific to individual particle time steps.
C=======================================================================
        INTEGER pactive,npactive,ntvector,itimestp,otimestp,inittbin,
     &          mintstep,upbin
        LOGICAL endstep
        REAL etol,stime,tsteppos

        COMMON/actcom/npactive,pactive(nbodsmax)
        COMMON/stepcom/itimestp(nbodsmax),otimestp(nbodsmax),inittbin,
     &                 mintstep,upbin,etol,stime,tsteppos,endstep,
     &                 ntvector

C=======================================================================
C   Definitions specific to vectorized tree construction, vectorized
C   tree walk, and vectorized tree search for nearest neighbors
C=======================================================================
        INTEGER pgroupb,subindex,bodlist,groups,templist,parent,nnear,
     &          asubp,celllist,groupbod,subpvect(nsubcell*ncells),
     &          ingroup,ngroups,npercell,nntot,nnmin,nnmax,nnavg,
     &          isubset,nworkvec
        REAL bottom,tempvect,workvect

        PARAMETER(ingroup=320,nworkvec=200000)

        EQUIVALENCE (subpvect(1),subp(nbods1,1))

        COMMON/nearcom/pgroupb(ncells+1),ngroups,groups(ncells),
     &                 groupbod(nbodsmax),npercell(nbods1:nbodcell),
     &                 nnear(nbodsmax),tempvect(nbodsmax),
     &                 bottom(1:nbodcell,ndim)
        COMMON/neighcom/nntot,nnmin,nnmax,nnavg
        COMMON/concom/celllist(ncells),parent(nbodsmax),asubp(nbodsmax),
     &                templist(nbodsmax),bodlist(nbodsmax),
     &                isubset(nbodsmax),subindex(nbodsmax)
        COMMON/workcom/workvect(nworkvec)



c-------------------------
c  hs 270307

        integer noutani,niraf
        real biraf
        common /hs/noutani,niraf,biraf      
