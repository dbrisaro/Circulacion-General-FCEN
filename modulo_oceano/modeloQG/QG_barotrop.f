      program QG_barotrop

c...............................................................................
c
c This program solves the barotropic vorticity equation in non-dimensional
c form using finite differences.
c The model has incorporated the "partial" slipping boundary conditions
c The wind has been set for a single gyre problem
c...............................................................................

c model parameters
c
c im: # of points in the x-direction
c jm: # of points in the y-direction
c ds: grid-step
c dt: time-step
c Ro: Rossby number. It is a measure of the non-linearity of the flow.
c eps: non-dimensional coefficient representing bottom friction
c Ah: nondimensional coefficient of horizontal Laplacian mixing
c Bh: nondimensional coefficient of horizontal bi-harmonic mixing
c gamma: coefficient of "intermediate slipping" used as boundary condition

c Arrays
c
c  p: stream function
c  psi: Laplacian of the stream function
c  curlt: wind stress curl
c  delpsi: delta 4th of the stream function. Used for horizontal mixing
c  Note: the subscript a,b,c denotes time steps t+1,t,t-1 respectively.

      double precision, allocatable, dimension (:,:) :: pa,pb,pc,psia,
     &  psib,psic,curlt,r,delpsi,delpsi4

      real kx, jpp, jxp, jpx, start, finish
      integer idx,BFP,MCF,GYR,HEM
      character*3 num_file
      logical :: chk_dir

      call cpu_time(start)

!      inquire(directory="out_tmp",exist=chk_dir)
!      print *,chk_dir
!      stop
      open(unit=17,file='./out_tmp/QG_diag.dat',
     *     form='formatted',status='unknown')
      open(unit=18,file='QG_param.dat',
     *     form='formatted',status='unknown')

      do i=1,7
        read(18,'(A)') dummy
      enddo
      read(18,'(A3,I8)') dump,im
      read(18,'(A3,I8)') dump,jm
      read(18,'(A3,F8.6)') dump,ds
      read(18,'(A3,F8.6)') dump,dt
      read(18,'(A3,F8.6)') dump,Ro
      read(18,'(A4,F8.6)') dump,eps
      read(18,'(A3,F8.6)') dump,Ah
      read(18,'(A3,F8.6)') dump,Bh
      read(18,'(A6,F8.6)') dump,gamma
      read(18,'(A4,I8)') dump,nst
      read(18,'(A5,I8)') dump,nend
      read(18,'(A5,I8)') dump,nplt
      read(18,'(A4,I2)') dump,MCF
      read(18,'(A6,I8)') dump,ncrit
      read(18,'(A6,F8.6)') dump,pcrit
      read(18,'(A4,I2)') dump,BFP
      read(18,'(A4,I2)') dump,GYR
      read(18,'(A4,I2)') dump,HEM
      close(18)

      allocate(pa(im,jm),pb(im,jm),pc(im,jm),psia(im,jm),psib(im,jm),
     &         psic(im,jm),curlt(im,jm),r(im,jm),delpsi(im,jm),
     &         delpsi4(im,jm))

c     general constants
      pi     = 3.1415169
      zero   = 0.0
      ele    = 2.*pi/(jm-1)
      kx     = pi/(im-1)
      dssqr  = ds**2
      dssqri = 1./dssqr
      dssqri4= dssqri/4.0
      epsdt  = eps*dt
      epsdti = 1./(1.+epsdt)
      rober  = 1.e-3
      nprts  = nprt
      nplts  = nplt
      c1     = epsdti
      c2     = BFP*(dt/ds)*epsdti
      c3     = 2.*dt*epsdti
      c4     = epsdt*epsdti
      c5     = 2.*dt*epsdti*Ro/3.0
      c6     = 2.*dt*Ah*epsdti
      c7     = 2.*dt*Bh*epsdti

      bc1 = -(1.-gamma*ds*0.5)/(1.+gamma*ds*0.5)
      bc2 =   1.-gamma*ds*0.5
      bc3 = - 1./(1. + gamma*ds*0.5)

      imm1=im-1
      imm2=im-2
      imm3=im-3
      imm4=im-4
      jmm1=jm-1
      jmm2=jm-2
      jmm3=jm-3
      jmm4=jm-4
 
c number of time steps to be saved
      xx= (nend-nst)/nplt
      nplots= aint(xx) + 1
      write(*,*) 'Time steps to be saved = ',nplots
      ksf=100
      kv=1000
      kp=0

c over-relaxation constants
c ncrit: # of steps allowed to do the relaxation in subroutine helm
c pcrit: criterium to stop the relaxation
c alfa: constant used for the sequential relaxation
c fxr:  constant used for the sequential relaxation

      const1= 1./(imm2-2)**2
      const2= 1./(jmm2-2.)**2
      alfa= 2.-4.442883*sqrt(const1+const2)
      fxr    = dssqr/4.0

c initialize all arrays

      do i=1,im
      do j=1,jm
        r      (i,j)= zero
        pa     (i,j)= zero
        pb     (i,j)= zero
        pc     (i,j)= zero
        psia   (i,j)= zero
        psib   (i,j)= zero
        psic   (i,j)= zero
        curlt  (i,j)= zero
        delpsi (i,j)= zero
        delpsi4(i,j)= zero
      enddo
      enddo

c wind stress curl

      if (GYR.eq.1) then ! single gyre

        do i=1,im
        do j=1,jm
          curlt(i,j)=-HEM*sin(ele*0.5*(j-1))
        enddo
        enddo

      elseif (GYR.eq.2) then ! double gyre

        do i=1,im
        do j=1,jm
          curlt(i,j)= -sin(ele*0.5*(j-1))*sin(kx*(i-1))
        enddo
        enddo

      else

        write(*,*) 'Wind Gyre type not defined'
        stop

      endif

      open(unit=20,file='./out_tmp/QG_wind_stress.dat',
     *     form='formatted',status='unknown') 
        do j=1,jm
          write(20,1000) (curlt(i,j),i=1,im)
        end do
      close(20)

c------------------------------------------------------------------------------
c
c                   T i m e    I n t e g r a t i o n
c
c------------------------------------------------------------------------------

      do ntime=nst,nend

!      write(*,*) 'step number = ',ntime

        do i=2,imm1
        do j=2,jmm1

c horizontal mixing

      delpsi(i,j)=(psic(i+1,j)+psic(i-1,j)+psic(i,j+1)+
     $                 psic(i,j-1)-4.*psic(i,j))*dssqri

      delpsi4(i,j)=(delpsi(i+1,j)+delpsi(i-1,j)+delpsi(i,j+1)+ 
     $                  delpsi(i,j-1)-4.*delpsi(i,j))*dssqri

c Arakawa's jacobian

      jpp= ((pb(i+1,j)-pb(i-1,j))*(psib(i,j+1)-psib(i,j-1))-
     $         (pb(i,j+1)-pb(i,j-1))*(psib(i+1,j)-psib(i-1,j)))
     $         *dssqri4

      jxp= (psib(i,j+1)*(pb(i+1,j+1)-pb(i-1,j+1))-psib(i,j-1)
     $          *(pb(i+1,j-1)-pb(i-1,j-1))-psib(i+1,j)*(pb(i+1,j+1)
     $          -pb(i+1,j-1))+psib(i-1,j)*(pb(i-1,j+1)-pb(i-1,j-1)))
     $          *dssqri4

      jpx=(pb(i+1,j)*(psib(i+1,j+1)-psib(i+1,j-1))-pb(i-1,j)
     $         *(psib(i-1,j+1)-psib(i-1,j-1))-pb(i,j+1)*(psib(i+1,j+1)
     $         -psib(i-1,j+1))+pb(i,j-1)*(psib(i+1,j-1)-psib(i-1,j-1)))
     $         *dssqri4

c update vorticity

      psia(i,j)= c1*psic(i,j) - c5*(jpp+jxp+jpx) - c2*(pb(i+1,j) -
     $               pb(i-1,j)) + c3*curlt(i,j) - c4*psic(i,j) + 
     $               c6*delpsi(i,j) - c7*delpsi4(i,j)

        enddo
      enddo

c update stream function

      nrelax = zero
100   nrelax = nrelax +1
      if(nrelax.gt.ncrit) goto 250

c solve the Laplacian for the stream function by over-relaxation

      do i=3,imm2
        do j=3,jmm2

         r(i,j)= ((pa(i-1,j)+pa(i+1,j)+pa(i,j-1)+pa(i,j+1)-4.*pa(i,j))*
     $             dssqri-psia(i,j))
         pa(i,j)= pa(i,j)+alfa*r(i,j)*fxr

        enddo
      enddo

c check convergence

      n1=0
      do i=3,imm2
        do j=3,jmm2

          rabs=sqrt(r(i,j)**2)
          if(rabs.gt.pcrit) n1=n1+1

        enddo
      enddo

      if(n1.ne.0) goto 100

      goto 300
250   write(*,*) 'Warning: The subroutine does not relax  ntime=',ntime
      goto 1122

300   continue


c set the boundary conditons on the stream function

      do j=2,jmm1
        pa(1,j )   = bc1*pa(3,j)
        pa(2,j)    = zero
        pa(im,j)   = bc1*pa(imm2,j)
        pa(imm1,j) = zero
      enddo

      do i=2,imm1
        pa(i,1)    = bc1*pa(i,3)
        pa(i,2)    = zero
        pa(i,jm  ) = bc1*pa(i,jmm2)
        pa(i,jmm1) = zero
      enddo

c diagonostic calculation of the vorticity on the walls

      do j=2,jmm1
        psia(2   ,j)= (pa(3 ,j) + pa(1   ,j))*dssqri
        psia(imm1,j)= (pa(im,j) + pa(imm2,j))*dssqri
      enddo

      do i=2,imm1
        psia(i,2)    = (pa(i,3 ) + pa(i,1   ))*dssqri
        psia(i,jmm1) = (pa(i,jm) + pa(i,jmm2))*dssqri
      enddo

c set the boundary conditons on the vorticity

      do j=2,jmm1
        psia(1,j)  = bc3*(bc2*psia(3,j)-4.*psia(2,j)+psia(2,j+1)+
     $               psia(2,j-1))
        psia(im,j) = bc3*(bc2*psia(imm2,j)-4.*psia(imm1,j)+
     $               psia(imm1,j+1)+psia(imm1,j-1))
      enddo

      do i=2,imm1
        psia(i,1)  = bc3*(bc2*psia(i,3)-4.*psia(i,2)+psia(i+1,2)+
     $               psia(i-1,   2))
        psia(i,jm) = bc3*(bc2*psia(i,jmm2)-4.*psia(i,jmm1)+
     $               psia(i+1,jmm1)+psia(i-1,jmm1))
      enddo

c time smoothing the stream function using a Robert's filter

      do i=1,im
        do j=1,jm
          pb  (i,j)= pb(i,j)+rober*(pa(i,j)-2.*pb(i,j)+pc(i,j))
          psib(i,j)= psib(i,j)+rober*(psia(i,j)-2.*psib(i,j)+psic(i,j))
        enddo
      enddo

c save data

c     kinetic energy (entire domain)

      tke=0.
      do i=2,imm1
        do j=2,jmm1
          gradpsi = ((pb(i+1,j)-pb(i-1,j))**2+(pb(i,j+1)
     $                 -pb(i,j-1))**2)*dssqri
          tke=tke+0.5*gradpsi
        enddo
      enddo

c     write local stream function and vorticity (domain center) and
c     total kinetic energy to file

      write(17,*) ntime,pb(im/2,jm/2),psib(im/2,jm/2),tke

      if(ntime.ge.nplts) then

        kp=kp+1
        ksf=ksf+kp
        kv=kv+kp
        write(num_file,'(I3)') kp
        if (kp.lt.10) then
          idx=3
        elseif (kp.ge.10.and.kp.lt.100) then
          idx=2
        else
          idx=1
        endif

        write(*,*) 'step number = ',ntime,' writing file = ',num_file

        if (idx.eq.3) then

          if(nplots.lt.100) then
            open(ksf,file='./out_tmp/psi0'//num_file(idx:3)//
     *      '.dat',form='formatted',status='unknown')
            open(kv,file='./out_tmp/vor0'//num_file(idx:3)//
     *      '.dat',form='formatted',status='unknown')
          else 
            open(ksf,file='./out_tmp/psi00'//num_file(idx:3)//
     *      '.dat',form='formatted',status='unknown')
            open(kv,file='./out_tmp/vor00'//num_file(idx:3)//
     *      '.dat',form='formatted',status='unknown')
          endif

        elseif (idx.eq.2) then

          if(nplots.lt.100) then
            open(ksf,file='./out_tmp/psi'//num_file(idx:3)//
     *      '.dat',form='formatted',status='unknown')
            open(kv,file='./out_tmp/vor'//num_file(idx:3)//
     *      '.dat',form='formatted',status='unknown')
          else
            open(ksf,file='./out_tmp/psi0'//num_file(idx:3)//
     *      '.dat',form='formatted',status='unknown')
            open(kv,file='./out_tmp/vor0'//num_file(idx:3)//
     *      '.dat',form='formatted',status='unknown')
          endif

        elseif (idx.eq.1) then

          open(ksf,file='./out_tmp/psi'//num_file(idx:3)//
     *    '.dat',form='formatted',status='unknown')
          open(kv,file='./out_tmp/vor'//num_file(idx:3)//
     *    '.dat',form='formatted',status='unknown')

        endif

        if (MCF.eq.0) then ! output in matrix form
          do j=1,jm
            write(ksf,1000) (pb(i,j),i=1,im)
          end do
        else ! output in columns
          do j=1,jm
            do i=1,im
              write(ksf,1001) i,j,pb(i,j)
            end do
          end do
        endif
        close(ksf)

        if (MCF.eq.0) then ! output in matrix form
          do j=1,jm
            write(kv,1000) (psib(i,j),i=1,im)
          end do
        else ! output in columns
          do j=1,jm
            do i=1,im
              write(kv,1001) i,j,psib(i,j)
            end do
          end do
        endif
        close(kv)

        nplts= nplts + nplt


      else

        write(*,*) 'step number = ',ntime

      endif

c update variables

      do i=1,im
        do j=1,jm
          pc  (i,j)= pb  (i,j)
          pb  (i,j)= pa  (i,j)
          psic(i,j)= psib(i,j)
          psib(i,j)= psia(i,j)
        enddo
      enddo

      enddo

 1000 format(1x,10000f10.3)
 1001 format(1x,2I4.0,f10.3)

 1122 close(16)
      close(17)

      call cpu_time(finish)
      write(*,*) 'Model run finished'
      print '("Elapsed time: ",f6.3," minutes")',(finish-start)/60

      end
