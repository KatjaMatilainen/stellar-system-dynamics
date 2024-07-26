c**************************************************************************

	subroutine outgrav

c**************************************************************************
	 INCLUDE 'treedefs.h'

         open(unit=28,file='outgrav')
         write(28,*) nbodies

         DO 70 p=1,nbodies
            write(28,*) mass(p)
 70      CONTINUE
         
         DO 80 p=1,nbodies
            write(28,*) pos(p,1),pos(p,2),pos(p,3)
 80      CONTINUE
         
         DO 90 p=1,nbodies
            write(28,*) vel(p,1),vel(p,2),vel(p,3)
 90      CONTINUE
         
         DO 230 p=1,nbodies
            write(28,*) acc(p,1),acc(p,2),acc(p,3)
 230     CONTINUE
         
         close(28)
         return
         end
         
