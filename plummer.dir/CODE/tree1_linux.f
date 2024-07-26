c     tree1.f
C***********************************************************************
C
C
                         FUNCTION ismax(n,x,inc)
C
C
C***********************************************************************
C
C
C     Function to locate index of maximum element of a real vector.
C
C
C=======================================================================

        INTEGER ismax,n,inc,i
        REAL x(1),xmax

        ismax=1
        xmax=x(1)

        DO 10 i=2,n,inc
           IF(x(i).GT.xmax) THEN
              ismax=i
              xmax=x(i)
           ENDIF
 10     CONTINUE

        RETURN 
        END

C***********************************************************************
C
C
                          FUNCTION ismin(n,x,inc)
C
C
C***********************************************************************
C
C
C     Function to locate index of minimum element of a real vector.
C
C
C=======================================================================

        INTEGER ismin,n,inc,i
        REAL x(1),xmin

        ismin=1
        xmin=x(1)

        DO 10 i=2,n,inc
           IF(x(i).LT.xmin) THEN
              ismin=i
              xmin=x(i)
           ENDIF
 10     CONTINUE

        RETURN 
        END

C***********************************************************************
C
C
                 FUNCTION isrchigt(n,iarray,inc,itarget)
C
C
C***********************************************************************
C
C
C     Function to return index of first element of iarray greater
C     than itarget, or n+1 if none is found.
C
C
C=======================================================================

        INTEGER isrchigt,n,iarray(1),inc,itarget,i

        isrchigt=n+1

        DO 10 i=1,n,inc
           IF(iarray(i).GT.itarget) THEN
              isrchigt=i
              RETURN
           ENDIF
 10     CONTINUE

        RETURN 
        END

C***********************************************************************
C
C
            SUBROUTINE wheneq(n,iarray,inc,itarget,index,nval)
C
C
C***********************************************************************
C
C
C     Subroutine to return locations of all elements of an integer
C     vector equal to itarget.
C
C
C=======================================================================

        INTEGER iarray(1),index(1),n,inc,itarget,nval,i

        nval=0

        DO 10 i=1,n,inc
           IF(iarray(i).EQ.itarget) THEN
              nval=nval+1
              index(nval)=i
           ENDIF
 10     CONTINUE

        RETURN 
        END

C***********************************************************************
C
C
            SUBROUTINE whenfgt(n,array,inc,target,index,nval)
C
C
C***********************************************************************
C
C
C     Subroutine to return locations of all elements of a real
C     vector greater than target.
C
C
C=======================================================================

        INTEGER index(1),n,inc,nval,i
        REAL array(1),target

        nval=0

        DO 10 i=1,n,inc
           IF(array(i).GT.target) THEN
              nval=nval+1
              index(nval)=i
           ENDIF
 10     CONTINUE

        RETURN 
        END

C***********************************************************************
C
C
            SUBROUTINE whenflt(n,array,inc,target,index,nval)
C
C
C***********************************************************************
C
C
C     Subroutine to return locations of all elements of a real
C     vector less than target.
C
C
C=======================================================================

        INTEGER index(1),n,inc,nval,i
        REAL array(1),target

        nval=0

        DO 10 i=1,n,inc
           IF(array(i).LT.target) THEN
              nval=nval+1
              index(nval)=i
           ENDIF
 10     CONTINUE

        RETURN 
        END

C***********************************************************************
C
C
            SUBROUTINE whenige(n,iarray,inc,itarget,index,nval)
C
C
C***********************************************************************
C
C
C     Subroutine to return locations of all elements of an integer
C     vector greater than or equal to itarget.
C
C
C=======================================================================

        INTEGER iarray(1),index(1),n,inc,itarget,nval,i

        nval=0

        DO 10 i=1,n,inc
           IF(iarray(i).GE.itarget) THEN
              nval=nval+1
              index(nval)=i
           ENDIF
 10     CONTINUE

        RETURN 
        END

C***********************************************************************
C
C
            SUBROUTINE whenigt(n,iarray,inc,itarget,index,nval)
C
C
C***********************************************************************
C
C
C     Subroutine to return locations of all elements of an integer
C     vector greater than itarget.
C
C
C=======================================================================

        INTEGER iarray(1),index(1),n,inc,itarget,nval,i

        nval=0

        DO 10 i=1,n,inc
           IF(iarray(i).GT.itarget) THEN
              nval=nval+1
              index(nval)=i
           ENDIF
 10     CONTINUE

        RETURN 
        END

C***********************************************************************
C
C
            SUBROUTINE whenile(n,iarray,inc,itarget,index,nval)
C
C
C***********************************************************************
C
C
C     Subroutine to return locations of all elements of an integer
C     vector less than or equal to itarget.
C
C
C=======================================================================

        INTEGER iarray(1),index(1),n,inc,itarget,nval,i

        nval=0

        DO 10 i=1,n,inc
           IF(iarray(i).LE.itarget) THEN
              nval=nval+1
              index(nval)=i
           ENDIF
 10     CONTINUE

        RETURN 
        END

C***********************************************************************
C
C
            SUBROUTINE whenilt(n,iarray,inc,itarget,index,nval)
C
C
C***********************************************************************
C
C
C     Subroutine to return locations of all elements of an integer
C     vector less than itarget.
C
C
C=======================================================================

        INTEGER iarray(1),index(1),n,inc,itarget,nval,i

        nval=0

        DO 10 i=1,n,inc
           IF(iarray(i).LT.itarget) THEN
              nval=nval+1
              index(nval)=i
           ENDIF
 10     CONTINUE

        RETURN 
        END

C***********************************************************************
C
C
            SUBROUTINE whenne(n,iarray,inc,itarget,index,nval)
C
C
C***********************************************************************
C
C
C     Subroutine to return locations of all elements of an integer
C     vector not equal to itarget.
C
C
C=======================================================================

        INTEGER iarray(1),index(1),n,inc,itarget,nval,i

        nval=0

        DO 10 i=1,n,inc
           IF(iarray(i).NE.itarget) THEN
              nval=nval+1
              index(nval)=i
           ENDIF
 10     CONTINUE

        RETURN 
        END
C***********************************************************************
C
C
                          SUBROUTINE ranset(sd)
C
C
C***********************************************************************
C
C
C     Dummy subroutine to initialize random numbers.
C
C
C=======================================================================

        REAL sd

        RETURN 
        END

C***********************************************************************
C
C
                           FUNCTION ranf(iran)
C
C
C***********************************************************************
C
C
C     Function to return random numbers.
C
C
C=======================================================================

        INTEGER iran
        REAL ranf,drand

        ranf=ran1(iran)


        RETURN 
        END

C***********************************************************************
C
C
                            FUNCTION second()
C
C
C***********************************************************************
C
C
C     Subroutine to return elapsed cpu time.
C
C
C=======================================================================

c        REAL etime,utime,stime,x,second
c        x=etime(utime,stime)
c        second=utime+stime     


          real tarray(2)
          double precision time
          real second

C     RETURNS THE ELAPSED TOTAL CPU TIME IN MSECS

          time=etime(tarray)
          second=time
   


        RETURN 
        END


      FUNCTION RAN1(IDUM)
      DIMENSION R(97)
      PARAMETER (M1=259200,IA1=7141,IC1=54773,RM1=3.8580247E-6)
      PARAMETER (M2=134456,IA2=8121,IC2=28411,RM2=7.4373773E-6)
      PARAMETER (M3=243000,IA3=4561,IC3=51349)


c      DATA IFF /0/
       save
c     IF (IDUM.LT.0.OR.IFF.EQ.0) THEN
      IF (IDUM.LT.0) THEN
       write(*,*) 'idum-seed',idum
        IFF=1
        IX1=MOD(IC1-IDUM,M1)
        IX1=MOD(IA1*IX1+IC1,M1)
        IX2=MOD(IX1,M2)
        IX1=MOD(IA1*IX1+IC1,M1)
        IX3=MOD(IX1,M3)
        DO 11 J=1,97
          IX1=MOD(IA1*IX1+IC1,M1)
          IX2=MOD(IA2*IX2+IC2,M2)
          R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
11      CONTINUE
c        write(*,*) r
      ENDIF
      IX1=MOD(IA1*IX1+IC1,M1)
      IX2=MOD(IA2*IX2+IC2,M2)
      IX3=MOD(IA3*IX3+IC3,M3)
      J=1+(97*IX3)/M3
c      IF(J.GT.97.OR.J.LT.1) PAUSE
      RAN1=R(J)
      R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1

c      write(99,*) 'rnd1: ',idum,ran1
      RETURN
      END


