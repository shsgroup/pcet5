c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*
c*      
c RCS ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rdiag(a,u,d,n,epsln)
c
c diagonalization of a real symmetric (n*n) matrix a by householder reduction
c
c INPUT:  a: real symmetric matrix stored in a linear array 
c            of dimension n*(n+1)/2, thus a(i,j) is  a(i+j*(j-1)/2) (j.ge.i)
c            a(ij) = <v(i)|a|v(j)> , v(i) denote basis vectors used to get a
c         u: n*n matrix with v(i) as column vectors
c         n: dimension of a
c         epslon: numerical threshhold, e.g. 1.d-14 for 64 bit words
c OUTPUT: u:  matrix of eigen vectors u(+)au = diagonal
c         d: eigenvalues
c 
c mainly taken from
c Press et al: "Numerical Recipes", Cambridge University Press, 1988
c     with appropriate modifications
c         as convenient for SCF procedures
c
c             R.A. Dec. 1990
c
c     with use of (level2) BLAS routines in all n**3 steps
c     Aug. 00  R A
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Revision : 2001/03/13 kakha
c  several scratch arraies replaced by the one allocatable scratch matrix 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit real*8 (a-h,o-z)
      dimension a(n*(n+1)/2),u(n,n),d(n)
      double precision e(n)
C     .. Scratch Array ..
      double precision, allocatable :: wspace(:,:)
c
c     quick return if possible
c
      if (n.le.0) then
        write (*,*) ' W A R N I N G: diagonalization called (rdiag)',
     * ' with dimension n = ',n
        return
      endif
      if (n.eq.1) then
        d(1)=a(1) 
        return
      endif
c
c      shift diagonals of a to make trace(a) = 0 (numerics)
c
      tr=0
      ii=0
      do i=1,n
         ii=ii+i
         tr=tr+a(ii)
c        d(i)=0.d0
      enddo
      shift=tr/dble(n)
      ii=0
      do i=1,n
         ii=ii+i
         a(ii)=a(ii)-shift
      enddo
c
c
c allocate Scratch Array wspace(:,:)
c
      if(.NOT.allocated(wspace)) allocate(wspace(n,1))
c ..
      call tredra (a,d,e,n,wspace,epsln)
c
c deallocate Scratch Array wspace(:,:)
c
      if(allocated(wspace)) deallocate(wspace)
c ..
c
c allocate Scratch Array wspace(:,:)
c
      iwork=40*n
      if(.NOT.allocated(wspace)) allocate(wspace(iwork,2))
      inwork=1+40*n
      call uta (a,u,n,wspace(1,1),wspace(1,2))
c
c deallocate Scratch Array wspace(:,:)
c
      if(allocated(wspace)) deallocate(wspace)
      nn=(n-1)*n/2+1
c
c allocate Scratch Array wspace(:,:)
c
      if(.NOT.allocated(wspace)) allocate(wspace(n,6))
      call tqli(d,e,u,n,epsln,wspace)
c
c deallocate Scratch Array wspace(:)
c
      if(allocated(wspace)) deallocate(wspace)
c
c     restore eigenvalues
c
      do 30 i=1,n
         d(i)=d(i)+shift
   30 continue
      return
      end
c*
c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*
c*
      subroutine tredra (a,d,e,n,dsc,epsln)
      implicit real*8 (a-h,o-z)
      dimension a(n*(n+1)/2),d(n),e(n),dsc(n)
c mainly taken from
c Press et al: "Numerical Recipes", Cambridge University Press, 1988
c     with appropriate modifications convenient for SCF procedures
c             R.A. Dec. 1990
c      Reduce a to tridiagonal form
c      Diagonals returned in vector d,
c      couplings in e 
c
c      modified 7-96 by R. A.
c      and again 8-00 by R. A. using BLAS routines in all n**3 steps
c      and again 3-01 by R. A. & Kakha :
c      now matrix-vector & matrix update operations are in one loop
c
      imi(i)=i*(i-1)/2
      a0=0.d0
      a1=1.d0
      a1m=-a1
      ah=0.5d0
c
c     special cases n.eq.2
c
      if (n.eq.2) then
        d(1)=a(1)
        e(1)=a(2)
        d(2)=a(3)
        return
      endif
c
c     get thr
c
      h=a0
      nn=imi(n+1)
      do ij=1,nn
         h=h+a(ij)**2
      enddo
      thr=epsln*sqrt(h/dble(nn))
      
      ii=imi(n)
      l=n-1
c	 
c    first form the vector d for the case i=n (if its norm .gt. thr)
c
            h=a0
            do k=1,l
               h=h+a(ii+k)**2
            enddo
            s=sqrt(h)
	    sn=s
         if (s.gt.thr) then
               f=a(ii+l)
               g=-sign(s,f)
               ei=g
               h=h-f*g
               a(ii+l)=f-g

               call dspmv('u',l,a1,a,a(ii+1),1,a0,d,1)
	   
               f=a0
               hi=a1/h
               do j=1,l
                  d(j)=d(j)*hi
                  f=f+d(j)*a(ii+j)
               enddo
               hh=f*hi*ah
               do j=1,l
                  d(j)=d(j)-hh*a(ii+j)
               enddo
         end if	       
c  ..
c
c     start reduction
c
      do i=n,4,-1
c
c     start reduction up to the case i=4 
c
         ii=imi(i)
         l=i-1

         if (sn.le.thr) then
            ei=a0
            do j=1,l
               a(ii+j)=a0
            enddo
         else
c
c     reduce the matrix a (if vector norm .gt. thr)
c
c     first do the case ij=l
c        if((d(l).ne.a0).or.(a(ii+l).ne.a0)) then
            term1 = -a(ii+l)
            term2 = -d(l)
	      ik = imi(l)
            do ji = 1, l 
               a(ik+ji)=a(ik+ji)+d(ji)*term1+a(ii+ji)*term2
            end do
c	 end if
c
c     end of case ij=l
c
         end if
c     form the vector u for the next i (if its norm .gt. thr)
c
            iin=imi(l)
            ln=l-1
         
            h=a0
            do k=1,ln
               h=h+a(iin+k)**2
            enddo
	    s=sn
            sn=sqrt(h)

         if (sn.gt.thr) then
            f=a(iin+ln)
            g=-sign(sn,f)
            ein=g
            h=h-f*g
            a(iin+ln)=f-g
	 end if   
c     ..	 
c 
c     next do the case ij=1,l-1 (if both vectors norm .gt. thr)
c
         if(s.gt.thr.AND.sn.gt.thr) then
            do ji = 1, ln
               dsc(ji)=a0
	    end do
	    ikk = 1
	    do ij = 1, l-1
               term1 = -a(ii+ij)
               term2 = -d(ij)
               term3 = a(iin+ij)
               term4 = a0
               ik = ikk
	       do ji = 1, ij-1 
                  a(ik+ji-1)=a(ik+ji-1)+d(ji)*term1+a(ii+ji)*term2
                  dsc(ji)=dsc(ji)+term3*a(ik+ji-1)
                  term4=term4+a(ik+ji-1)*a(iin+ji)		 
               end do
               a(ik+ji-1)=a(ik+ji-1)+d(ij)*term1+a(ii+ij)*term2
               dsc(ij)=dsc(ij)+term3*a(ikk+ij-1)+term4
               ikk = ikk + ij
            end do
            do ji = 1, ln
               d(ji)=dsc(ji)
            end do
	 else if(sn.gt.thr) then
            do ji = 1, ln
               dsc(ji)=a0
	    end do
            call dspmv('u',ln,a1,a,a(iin+1),1,a0,dsc,1)	    
            do ji = 1, ln
               d(ji)=dsc(ji)
            end do
         else if(s.gt.thr) then
	    ikk = 1
	    do ij = 1, l-1
               term1 = -a(ii+ij)
               term2 = -d(ij)
               ik = ikk
	       do ji = 1, ij-1 
                  a(ik+ji-1)=a(ik+ji-1)+d(ji)*term1+a(ii+ji)*term2
               end do
               a(ik+ji-1)=a(ik+ji-1)+d(ij)*term1+a(ii+ij)*term2
               ikk = ikk + ij
            end do
	 end if
c	    
c     end of case ij=1,l-1
c
c     form the vector d for the next i (if its norm .gt. thr)
c
         if (sn.gt.thr) then
	    f=a0
            hin=a1/h
            do j=1,ln
               d(j)=d(j)*hin
               f=f+d(j)*a(iin+j)
            enddo
            hh=f*hin*ah
            do j=1,ln
               d(j)=d(j)-hh*a(iin+j)
            enddo
         end if   
c     ..	    

         if(s.gt.thr) then
            hinv=sqrt(hi)
            do j=1,l
               a(ii+j)=a(ii+j)*hinv
            enddo
         end if   
c
c     end of reduction up to the case i=4
c
         d(i)=a(ii+i)
         e(l)=ei

         ei=ein
         hi=hin
      enddo
c
c     reduce the matrix a for the last case i=3
c
         ii=imi(3)
         l=2

         if (sn.le.thr) then
            ei=a0
            do j=1,l
               a(ii+j)=a0
            enddo
         else
            call dspr2('u',l,a1m,d,1,a(ii+1),1,a)
	    
            hinv=sqrt(hi)
            do j=1,l
               a(ii+j)=a(ii+j)*hinv
            enddo
	 end if   

         d(3)=a(ii+3)
         e(l)=ei	    
c
c     end of reduction for the case i=3
c
c     end of reduction
c
      e(1)=a(2)
      d(1)=a(1)
      d(2)=a(3)
      return
      end

c*
c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*
c*
      subroutine uta(a,u,n,d,as)
      implicit real*8 (a-h,o-z)
      dimension a(n*(n+1)/2),u(n,n),d(n,40),as(n,40)
      imi(i)=i*(i-1)/2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Purpose
c     =======
c     matrix multiply  u*v
c     where v is an orthogonal matrix defined by a,
c     which corresponds to successive reflections
c     as usual in Householder procedures.
c     A straight fortran code would be:
c     do i=n,3,-1
c        ii=imi(i)
c        l=i-1
c        call dgemv ('n',n,l,a1,u,n,a(ii+1),1,a0,d,1)
c        call dger(n,l,a1m,d,1,a(ii+1),1,u,n)
c     enddo
c
c     Arguments
c     =========
c
c     a  (input) double precision array, dimension (n*(n+1)/2).
c        On entry, corresponds to successive reflections,packed   
c        columnwise in a linear array, as returned by tredra.
c        
c     u  (input/output) double precision array, dimension(n,n).
c        On entry, an orthogonal matrix from previuse iteration.
c        On exit, overwritten by matrix-matrix maltiple.
c
c     n  (input) The order of the symmetric tridiagonal matrix.
c      
c     d  (workspace) double precision array, dimension(n,40)
c
c     as (workspace) double precision array,dimension(n,40)
c
c     Dependences
c     ===========
c     uses level 3 Blas routine dgemm
c
c     Further Details
c     ===============
c     uses the optimized block size : iblk=40
c     
c     Mar. 2001 R. A.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      a0=0.d0
      a1=1.d0
      a1m=-a1
      iblk=40
      imin=44
      l=n-1+iblk
      do i=n,imin,-iblk
        l=i-1
c        copy strips of a on matrix as
        kend=i-iblk
        do ib=1,iblk
          iib=imi(kend+1)
          do k=1,kend
            as(k,ib)=a(iib+k)
          enddo
          do k=kend+1,l
            as(k,ib)=a0
          enddo
          kend=kend+1
        enddo
c       end copy

        call dgemm('n','n',n,iblk,l,a1,u,n,as,n,a0,d,n)

c       correct as
        ke=l-1
        do j=iblk-1,1,-1
          scal=a0
          do k=1,ke
            scal=scal+as(k,j)*as(k,j+1)
          enddo
          do jj=j+1,iblk-1
            scaln=a0
            do k=1,ke
              as(k,jj)=as(k,jj)-scal*as(k,j)
              scaln=scaln+as(k,j)*as(k,jj+1)
            enddo
          scal=scaln
          enddo
          do k=1,ke
            as(k,iblk)=as(k,iblk)-scal*as(k,j)
          enddo
          ke=ke-1
        enddo
c       end correct as

        call dgemm('n','t',n,l,iblk,a1m,d,n,as,n,a1,u,n)
      enddo

c     process remaining cases
      irem=l-iblk+1
      do i=irem,3,-1
         ii=imi(i)
         l=i-1
         call dgemv ('n',n,l,a1,u,n,a(ii+1),1,a0,d,1)
         call dger(n,l,a1m,d,1,a(ii+1),1,u,n)
      enddo
      return
      end
c
c*
c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*
c*
      subroutine trtbls(n,m,x,b,c,v)
      implicit real*8(a-h,o-z)

c =====================================================================
c
c     x = transpose(c) * b * c
c
c =====================================================================
c    using BLAS routines, which necessitates scratch arrays
c    scr1(n,100),scr2(n,100) holding columns of symmetric matrices
c    blocking is fixed at 100
c

      dimension x(m*(m+1)/2),b(n*(n+1)/2),c(n,m)
     *  ,scr1(m,100),scr2(m,100)

      md=100
      a1=1.d0
      a0=0.d0
      do i=1,md
        do j=1,n
          scr1(j,i)=a0
          scr2(j,i)=a0
        enddo
      enddo
      is=1
      ie=min(m,md)
      nblk=(m+md-1)/md
      do iblk=1,nblk
        nrws=ie-is+1
        call smdgm (scr2,b,c(1,is),scr1,n,md,nrws)
        call dgemm ('t','n',ie,nrws,n,a1,c,n,scr2,n,a0,scr1,m)
        call putstrs (scr1,x,is,ie,m,md)
        is=ie+1
        ie=min(ie+md,m)
      enddo
      return
      end
c*
c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*
c*
      subroutine getstrp (as,a,is,ie,n,nbmax)
      implicit real*8(a-h,o-z)
      dimension a(n*(n+1)/2),as(n,nbmax)
c
c     Get a strip of columns is..ie out of the symmetric matrix a
c     of dimension n, a is in (upper) packed form.
c     This is put on the matrix as in the first (ie-is+1) columns
c     where a has dimensions (n,nbmax)
c     Aug 00 R A
c
      if (is.lt.1.or.ie.lt.is.or.ie.gt.n) stop 'getstrp: is or ie'
      jj=is*(is-1)/2
      ioff=is-1
      do j=is,ie
        do i=1,j
          as(i,j-ioff)=a(jj+i)
        enddo
        jj=jj+j
      enddo
      do j=is,ie
        do i=j+1,ie
          as(i,j-ioff)=as(j,i-ioff)
        enddo
      enddo
      ii=ie*(ie+1)/2
      do i=ie+1,n
        do j=is,ie
          as(i,j-ioff)=a(ii+j)
        enddo
        ii=ii+i
      enddo
      return
      end
c*
c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*
c*
      subroutine putstrs (as,a,is,ie,n,nbmax)
      implicit real*8(a-h,o-z)
      dimension a(n*(n+1)/2),as(n,nbmax)
c
c     Purt a strip of columns is..ie of the matrix as
c     of dimension (n,*) onto the corresponding columns
c     of the symmetric matrix a, which is in (upper) packed form
c     Aug 00 R A
c
      if (is.lt.1.or.ie.lt.is.or.ie.gt.n) stop 'putstr: is or ie'
      jj=is*(is-1)/2
      ioff=is-1
      do j=is,ie
        do i=1,j
          a(jj+i)=as(i,j-ioff)
        enddo
        jj=jj+j
      enddo
      return
      end
c*
c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*
c*
      subroutine smdgm (a,s,b,scr,n,md,nrow)
      implicit real*8 (a-h,o-z)
      dimension a(n,md),b(n,md),scr(n,md),s(n*(n+1)/2)
c
c     matrix multiply  a = s * b
c     s is symmetric packed of dimension n
c     b and a are (n*nrow) matrices
c     scr  is a scratch array (n*md).
c     Using BLAS routines which necessitates to get and put
c     columns of (anti) symmetric matrices on a normal matrix
c     Aug 00 R A
c
      a1=1.d0
      a0=0.d0
      is=1
      ie=min(n,md)
      nblk=(n+md-1)/md
      do iblk=1,nblk
        nrws=ie-is+1
        call getstrp (scr,s,is,ie,n,nrws)
        call dgemm ('t','n',nrws,nrow,n,a1,scr,n,b,n,a0,a(is,1),n)
        is=ie+1
        ie=ie+md
        ie=min(ie,n)
      enddo
      return
      end
c*
c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*
c*
c RCS ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c $Id: rdiag.f,v 5.1 2011-04-13 23:49:48 souda Exp $
c $Log: not supported by cvs2svn $
c Revision 1.4  2001/02/19 13:00:20  peter
c Changes to get NOs converged in big molecules
c
c Revision 1.3  2000/10/20 12:55:57  haettig
c replaced stop commandos by a call to quit.
c
c Revision 1.2  2000/08/31 16:53:00  uwe
c diagonalization now using BLAS Routines, it's a two step procedure: in the
c first sweep only eigenvalues are computed, in the 2nd sweep they are used
c to get the vectors .........                   by R.A., 31.08.2000
c
c Revision 1.1  1997/02/18 12:26:58  marco
c *** empty log message ***
c
c Revision 1.2  1996/08/28 14:45:00  karin
c new tqli: R.A. 7.1996
c
c Revision 1.1  1992/09/10 15:00:24  tomjones
c Initial revision
c
c RCS ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine tqli(d,e,u,n,epsln,work)
      implicit real*8 (a-h,o-z)
      dimension u(n,n),d(n),e(n)
      dimension work(n,6)
      logical unprop
c
c mainly taken from
c Press et al: "Numerical Recipes", Cambridge University Press, 1988
c   with appropriate modifications convenient for SCF procedures
c             R.A. Dec. 1990
c     With use of (level2) BLAS routines in all n**3 steps.
c     A two-sweep procedure is used:
c         in the first sweep only eigenvalues are computed,
c         these are used in the 2nd sweep to get also vectors,
c         this requires only about a single iteration to converge
c         since the eigenvalue is already exact .
c     Aug. 00  R A
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      unprop=.false.
cdebug      istat1=0
cdebug      istat2=0
      a0=0.d0
      a1=1.d0
      if (n.le.1) return
c      get thr
      thr=abs(d(1))
      do i=2,n
         thr=max(thr,abs(d(i)),abs(e(i-1)))
      enddo
      thr0=epsln*thr
      thr1=max(2.d-15,epsln)*thr
      e(n)=a0
c     save d and e
      do i=1,n
        work(i,1)=d(i)
        work(i,2)=e(i)      
      enddo
      do iswp =1,2
      nprim=0
      ngivrt=0
      irot=0
      do l=1,n
         thr=thr0
         iter=0
    1    do m=l,n-1
            if (abs(e(m)).le.thr) go to 2
         enddo
         m=n
    2    if (m.ne.l) then

            if (iter.gt.15) then
              thr=thr1
              if (iter.eq.25) thr=thr+thr
            endif
            if (iter.ge.300) 
     &        call quit( ' too many iterations in tqli (rdiag)')

            ngivrt=1+ngivrt
            iter=iter+1
cdebug            if (iswp.eq.2 .and. iter.gt.2.and.(.not.unprop)) then
cdebug              write(*,*) ' tqli: problem with ordering'
cdebug            endif
            if (iswp.eq.2 .and. iter.gt.2) unprop=.true.
            if ((iswp.eq.1).or.unprop) then
              g=(d(l+1)-d(l))/(e(l)+e(l))
              r=sqrt(g*g+a1)
              g=d(m)-d(l)+e(l)/(g+sign(r,g))
            else
	      g=d(m)-work(l,2)
            endif
            s=a1
            c=a1
            p=a0
            do i=m-1,l,-1
               nprim=1+nprim
               f=s*e(i)
               b=c*e(i)
               if (abs(f).ge.abs(g)) then
                  c=g/f
                  r=sqrt(c*c+a1)
                  e(i+1)=f*r
                  s=a1/r
                  c=c*s
               else
                  s=f/g
                  r=sqrt(s*s+a1)
                  e(i+1)=g*r
                  c=a1/r
                  s=s*c
               endif
               g=d(i+1)-p
               r=(d(i)-g)*s+2*c*b
               p=s*r
               d(i+1)=g+p
               g=c*r-b
               work(i,3)=c
               work(i,4)=s	       
            enddo
            if (iswp.eq.2) then
              irot=irot+1
              if (irot.eq.1) then
                do k=l,m-1
                  work(k,5)=work(k,3)
                  work(k,6)=work(k,4)		
                enddo
                mo=m
                lo=l
              else
                if (mo.eq.m.and.lo+1.ge.l) then
		  call rotall2(n,u,work(1,3),n,m,l,lo)
cdebug                  istat2=istat2+1
                else
cdebug                  istat1=istat1+1
                  call rotall(n,u,work(1,5),n,mo,lo)
                  call rotall(n,u,work(1,3),n,m,l)
                endif
                irot=0
              endif
            endif
            d(l)=d(l)-p
            e(l)=g
            e(m)=a0
            go to 1
         endif
      enddo
      if (iswp.eq.1) then
        do i=1,n
          tmp=work(i,2)
          work(i,2)=d(i)
          d(i)=work(i,1)	
          e(i)=tmp
        enddo
      endif
c      thr=thr*2.d0
      enddo
      if (irot.eq.1) call rotall(n,u,work(1,5),n,mo,lo)      
cdebug      if (irot.eq.1) istat1=istat1+1
      return
      end
c
      subroutine rotall (n,u,space,ilds,m,l)
      implicit real*8 (a-h,o-z)
      dimension u(n,n),space(ilds,*)
      
      if ((m-1).lt.l) return
      mmo=mod(m-l,2)
      if (mmo.eq.1) then
        i=m-1
        do j=1,n
          r0=u(j,i+1)
          r1=u(j,i)
          u(j,i+1)=space(i,2)*r1+space(i,1)*r0
          u(j,i)=space(i,1)*r1-space(i,2)*r0	  
        enddo
      endif
c
      do i=m-1-mmo,l,-2
        do j=1,n
          r0=u(j,i+1)
          s0=u(j,i)
          s1=u(j,i-1)
          r1=space(i,1)*s0-space(i,2)*r0
          u(j,i+1)=space(i,2)*s0+space(i,1)*r0
          u(j,i)=space(i-1,2)*s1+space(i-1,1)*r1
          u(j,i-1)=space(i-1,1)*s1-space(i-1,2)*r1	  
        enddo
      enddo
      return
      end
c
      subroutine rotall2 (n,u,space,ilds,m,l,lo)
      implicit real*8 (a-h,o-z)
      dimension u(n,n),space(ilds,*)
      
      i=m-1
      do j=1,n
        r0=u(j,i+1)
        r1=u(j,i)
        u(j,i+1)=space(i,4)*r1+space(i,3)*r0
        u(j,i)=space(i,3)*r1-space(i,4)*r0
      enddo
c
      do i=m-1,lo+1,-1
        do j=1,n
          r0=u(j,i+1)
          r1=u(j,i)
          r2=u(j,i-1)
          s2=space(i-1,4)*r2+space(i-1,3)*r1
          u(j,i-1)=space(i-1,3)*r2-space(i-1,4)*r1
          u(j,i)=space(i,1)*s2-space(i,2)*r0
          u(j,i+1)=space(i,2)*s2+space(i,1)*r0	  
        enddo
      enddo
      if (lo.eq.l) then
        i=l
        do j=1,n
          r0=u(j,i+1)
          r1=u(j,i)
          u(j,i+1)=space(i,2)*r1+space(i,1)*r0
          u(j,i)=space(i,1)*r1-space(i,2)*r0
        enddo
      endif
      return
      end
c*
c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*
c*
*======================================================================*
      subroutine quit(error_message)
*----------------------------------------------------------------------*
*
*  Purpose: abort program without taking care of MPI environment
*
*----------------------------------------------------------------------*
      implicit none
 
      common /thisis/ thisis
      character*80    thisis
      character*(*) error_message
      integer ll, lstr

c      ll = lstr(thisis,80)
       ll = 80

c print error message into output file
      write(6,'(/1x,a)') error_message
      write(6,'(1x,a)') thisis(1:ll)//' ended abnormally'

c flush and close output file
      close(6)

c notify poor user that something has happened
c      call cerr(thisis(1:ll)//' ended abnormally'//char(0))
      write(6,'(1x,a)') thisis(1:ll)//' ended abnormally'
      stop 

      end
c*
c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*
c*

