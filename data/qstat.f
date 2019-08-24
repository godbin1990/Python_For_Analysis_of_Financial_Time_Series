c
c     Computes portmanteau statistic for multivariate time series
c     regression for a ik-dimensional process.
c
c    The program has been checked, but some errors may exist. 
c    It is intended for use in my courses only. Ruey S. Tsay, December 1990.
c
      implicit double precision(a-h,o-z)
      parameter(mxn=100000, id=10, mlag = 120)

      dimension e(mxn,id), cov(id,id), corr(id,id), 
     *     trcorr(id,id), sigma(id,id), esq(mxn,id)
      dimension siginv(id,id) 

c     Read in data
c
      call setup(e,mxn,nob,id,ik,nlag,npar)
c
c  remove sample means
c
      do i=1, ik
       ave  = 0.0d0
       do  j=1, nob
        ave = ave + e(j,i)
       enddo
       ave = ave/dfloat(nob)
       do j=1, nob
        e(j,i) = e(j,i) - ave
       enddo
      enddo
c
      write(6,41)
 41   format(1x,'Compute Q-stat for the squared series? (1=y): ',$)
      read(5,*) isq
c
      if(isq .ne. 1)go to 45
c------------- square series
      do i=1, ik
       do it=1, nob
        esq(it,i) = e(it,i)*e(it,i)
       enddo
      enddo
c
c  remove sample means
c
      do i=1, ik
       ave  = 0.0d0
       do  j=1, nob
        ave = ave + esq(j,i)
       enddo
       ave = ave/dfloat(nob)
       do j=1, nob
        esq(j,i) = esq(j,i) - ave
       enddo
      enddo
c
 45   continue
      print*,'Output file: fort.36'

      tt = dlog(dfloat(nob))
      idlag=ifix(tt)

      write(36,*)'Suggested # of Lags = ',idlag
c      print*,'Idea # of lags = ',idlag
c
c     Compute sample covariance matrix

      call covar(e,sigma,mxn,nob,id,ik,0)
c
      call mtinv(sigma,siginv,id,ik)
c
      q = 0.0d0
      qcross = 0.0d0

      iend = nlag
      if(iend .gt. mlag)iend= mlag
c
      write(36,38)
 38   format(2x,'m',7x,'Orig-Q',7x,'LM-Q',7x,'LB-Q',5x,'d.f.')
      do it = 1, iend
        call covar(e,cov,mxn,nob,id,ik,it)
        call correlation(cov,corr,trcorr,siginv,id,ik)
        call mult(corr,trcorr,trace,trcross,mxn,nob,id,ik,it)
        q = q + trace
        qcross = qcross + trcross
        qstar = q + dfloat((ik**2)*it*(it+1))/(2.0d0*dfloat(nob))
        idf = it*ik*ik-npar
        if(idf .lt. 0)idf = 0
c
        write(36,37) it, q, qstar, qcross, idf
      enddo
   37   format(1x,i3,1x,3f12.4,i5)
cc****
c*** squared series
c****
      if(isq .ne. 1)stop
       write(36,42)
 42    format(1x,'*** For squared series *** ')
c
c     Compute sample covariance matrix of the squared series

      call covar(esq,sigma,mxn,nob,id,ik,0)
c
      call mtinv(sigma,siginv,id,ik)
c
      q = 0.0d0
      qcross = 0.0d0
c
      write(36,38)
      do it = 1, iend
        call covar(esq,cov,mxn,nob,id,ik,it)
        call correlation(cov,corr,trcorr,siginv,id,ik)
        call mult(corr,trcorr,trace,trcross,mxn,nob,id,ik,it)
        q = q + trace
        qcross = qcross + trcross
        qstar = q + dfloat((ik**2)*it*(it+1))/(2.0d0*dfloat(nob))
        idf = it*ik*ik-npar
        if(idf .le. 0)idf = 0
c
        write(36,37) it, q, qstar, qcross, idf
      enddo
c
      stop
      end

c------------------------------------------------------------------
      subroutine covar(resid,c,n,nob,id,ik,lag)
      implicit double precision(a-h,o-z)
c
c     Computes covariance matrices

      dimension resid(n,id), c(id,id)

      do i = 1,ik
        do j = 1,ik
          c(i,j) = 0.0d0
        enddo
      enddo

      do i = 1,ik
        do j = 1,ik
          accu=0.0d0
          do k = lag + 1, nob
            accu = accu + resid(k,i)*resid(k-lag,j)
          enddo
          c(i,j) = accu
        enddo
      enddo      

c     Divide matrix by number of observations, n

      do i = 1,ik
        do j = 1,ik
          c(i,j) = c(i,j)/dfloat(nob)
        enddo
      enddo

      return
      end

c------------------------------------------------------------------

      subroutine correlation(c,r,rt,sinv,id,ik)
      implicit double precision(a-h,o-z)
c
c     Computes correlation matrices

      dimension c(id,id), r(id,id), rt(id,id), sinv(id,id)

c     Compute r = Rk

      do i = 1,ik
        do j = 1,ik
          r(i,j) = 0.0d0
        enddo
      enddo

      do i = 1,ik
        do j = 1,ik
          do k = 1,ik
c            r(i,j) = r(i,j) + sinv(i,k)*c(k,j)
            r(i,j) = r(i,j) + c(i,k)*sinv(k,j)
          enddo
        enddo
      enddo

c     Compute rt = R-k

      do i = 1,ik
        do j = 1,ik
          rt(i,j) = 0.0d0
        enddo
      enddo

      do i = 1,ik
        do j = 1,ik
          do k = 1,ik 
c            rt(i,j) = rt(i,j) + sinv(i,k)*c(j,k)
            rt(i,j) = rt(i,j) + c(k,i)*sinv(k,j)
          enddo
        enddo
      enddo

      return
      end

c------------------------------------------------------------------

      subroutine mult(r,rt,sum,sumcr,n,nob,id,ik,lag)
c
      implicit double precision(a-h,o-z)
c
c     Multiplies r and rt, computes trace, then multiplies by n.
c
      parameter(m=10)
      dimension r(id,id), rt(id,id), temp(m,m)

c     Multiply r and rt:

      do i = 1,ik
        do j = 1,ik
          temp(i,j) = 0.0d0
        enddo
      enddo

      do i = 1,ik
        do j = 1,ik
          do k = 1,ik
c            temp(i,j) = temp(i,j) + r(i,k)*rt(k,j)
            temp(i,j) = temp(i,j) + rt(i,k)*r(k,j)
          enddo
        enddo
      enddo

c     Compute trace

      sum = 0.0d0

      do i = 1,ik
        sum = sum + temp(i,i)
      enddo

      sumcr = (dfloat(nob*nob))*(sum/(dfloat(nob-lag)))

      sum = sum*dfloat(nob)

      return
      end

c---------------------------------------------------------------

      subroutine setup(resid,n,nob,id,k,nlag,npar)
      implicit double precision(a-h,o-z)

c     Read in residuals

      dimension resid(n,id), jdx(10)
      character*40  fname
c
      write(6,1)
    1 format(1x,'residual (data) file_name: ',$)
      read(5,2) fname
    2 format(a40)

      write(6,3)
    3 format(1x,'#(columns) in the data file: ',$)
      read(5,*) nc
      write(6,34)
 34   format(1x,'#(parameters) if a model was fit: ',$)
      read(5,*) npar
      write(6,35)
 35   format(1x,'maximum number of CCM mtx: ',$)
      read(5,*) nlag
      if(nlag.le.0)nalg = 12
      if(npar.le.0)npar = 0
c
      open(21,file=fname,status='old')

      nob = 0
      do 10 i = 1,n
        read (21,*,end=11) (resid(i,j), j=1,nc)
         nob = nob+1
 10     continue
 11     close(21)
c
       print*,'nob = ', nob
c
      write(6,4)
 4    format(1x,'The number of series: ',$)
      read(5,*) k

      if(k.eq.nc)return
      write(6,5) k
 5    format(1x,'Input',i3,' locators of residual series: ',$)
      read(5,*) (jdx(i),i=1,k)
c
      do 15 i=1, nob
       do 13 j=1, k
 13       resid(i,j) = resid(i,jdx(j))
 15       continue
c
      return
      end
C*******************************************************************
      SUBROUTINE MTINV(A,DA,IDM,IDIM)
C***********************************************************************
C     MTINV COMPUTES THE INVERSE MATIRX OF A
C**********************************************************************
C**** 
      implicit double precision(a-h,o-z)
      dimension A(IDM,IDM),DA(IDM,IDM)
C****
      DETA=1.0D0
      IF(IDIM .EQ. 1) GO TO 600                                    
C****                                                                         
      DO 100 I=1,IDIM                                                         
      DO 100 J=1,IDIM                                                         
  100 DA(I,J)=A(I,J)                                                          
C****                                                                         
  120 DO 500 I=1,IDIM                                                         
      PIVOT=DA(I,I)                                                           
      DETA=DETA*PIVOT                                                         
C****                                                                         
C**** DIVIDE PIVOT ROW BY PIVOT ELEMENT.                                      
C****                                                                         
      DA(I,I)=1.0D0
      DPIVOT=PIVOT+1.0D-20                                              
      DPIVOT=DA(I,I)/DPIVOT                                             
      PIVOT=DPIVOT                                                      
      DO 200 J=1,IDIM                                                         
  200 DA(I,J)=DA(I,J)*PIVOT                                                   
C****                                                                         
C**** REDUCE NON-PIVOT ROWS                                                   
C****                                                                         
  210 DO 500 II=1,IDIM                                                        
      IF(II.EQ.I) GO TO 500                                                   
      T=DA(II,I)                                                              
      DA(II,I)=0.0D0
      DO 300 J=1,IDIM                                                         
  300 DA(II,J)=DA(II,J)-DA(I,J)*T                                             
C****                                                                         
  500 CONTINUE                                                          
      RETURN  
  600 DA(1,1)=DETA/A(1,1)
      RETURN
      END                                                                     
C****                                                                         
C**** END OF 'MTINV'                                                          
C****                                                                         
C***********************************************************************      
