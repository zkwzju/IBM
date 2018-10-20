      subroutine main_run
      use common_m
      use init_m
      REAL*8 espz(nz) 
      REAL*8 espx(nxh,nz),dspx(nxh,nz),espy(nyh,nz),dspy(nyh,nz)
      REAL*8 fmean(nz),xrms(nz),yrms(nz),pmean(nz),prms(nz)
      REAL*8 tfmean,txrms,tyrms,tpmean,tprms !,igmod
      REAL*8 auxdepo(nx0,ny0),auxeros(nx0,ny0)

      call initialization
c     loop over the time steps
      do itime = itfirst,iend

c     write the time step number to logfile
        if (mod((itime-istart),oi_cfl).eq.0) then
          write(26,*)' '
          write(26,*)' itime =  ',itime
          write(26,*)' '
        endif
c
c     perform one complete time step
        call advance(0)
        rtime  = rtime+dt
c
c     compute mean values
!>>>
!     call make_means
!<<<
c NOTE: adcp data for JR is output in subroutine 
c       output_adcpdata that is called in make_means above
c
c     output means
        if(FLAVER)then
c     keep these two "if" like this because oi_mean=0 -> FLAVER=.false.
           if(mod((itime-istart),oi_mean).eq.0)call output_means
        endif
c
c     compute deposit, must be done every time step        
        if(.false..and.FLSTLE)then
           call dcopy(nx0y0,tt(1,1,nz),1,auxdepo,1) !copy bottom concentration
           call fft_f2r_2d_new_slice(auxdepo) !bring it to real space
           call comp_erosion(auxeros) !erosion returns in real space
           do j=1,ny
             do i=1,nx
                erosion(i,j)= erosion(i,j)-ws*auxeros(i,j)*dt
                deposit(i,j)= deposit(i,j)+ws*auxdepo(i,j)*dt
             enddo
           enddo
c
        endif
c
c     produce outputs in Fourier space
c     output spectra
        if((mod(itime-istart,oi_spec).eq.0.or.
     &      mod(itime-istart,oi_mean).eq.0)) then
c u 
           call make_spec_new(u,espx,dspx,espy,dspy,xrms,
     &          txrms,yrms,tyrms)
c     if(mod(itime,oi_mean).eq.0)call make_meanf(u,fmean,tfmean)
           if(mod(itime-istart,oi_spec).eq.0)
     &          call spec_out(espx,dspx,espy,dspy,'dspu',itime)
c v
           call make_spec_new(v,espx,dspx,espy,dspy,xrms,
     &          txrms,yrms,tyrms)
c     if (mod(itime,oi_mean).eq.0)call make_meanf(v,fmean,tfmean)
           if (mod(itime-istart,oi_spec).eq.0)
     &          call spec_out(espx,dspx,espy,dspy,'dspv',itime)
c w
           call make_spec_new(w,espx,dspx,espy,dspy,xrms,
     &          txrms,yrms,tyrms)
c     if (mod(itime,oi_mean).eq.0)call make_meanf(w,fmean,tfmean)
           if (mod(itime-istart,oi_spec).eq.0)
     &          call spec_out(espx,dspx,espy,dspy,'dspw',itime)
c tt
           if (FLTHRM) then
              call make_spec_new(tt,espx,dspx,espy,dspy,xrms,
     &             txrms,yrms,tyrms)
c     if (mod(itime,oi_mean).eq.0)call make_meanf(tt,fmean,tfmean)
              if (mod(itime-istart,oi_spec).eq.0)
     &             call spec_out(espx,dspx,espy,dspy,'dspT',itime)

           endif
        endif
c
c     produce outputs in real space
        if( mod((itime-istart),oi_chan).eq.0.or.
     &      mod((itime-istart),oi_cfl ).eq.0.or.
     &      mod((itime-istart),oi_mean).eq.0.or.
     &      mod((itime-istart),oi_spec).eq.0.or.
     &      mod((itime-istart),oi_gbal).eq.0.or.
     &      mod((itime-istart),oi_1d  ).eq.0)then
c
c     transform to real space
           call fft_f2r_2d_new(u)
           call fft_f2r_2d_new(v)
           call fft_f2r_2d_new(w)
           if (FLTHRM)call fft_f2r_2d_new(tt)
c NOTE: old vel are not fft, they are always in Fourier space
c
c     output a restart file
           if(mod((itime-istart),oi_chan).eq.0) call output
c
c     output pressure (press is in Fourier space in var p allways)
           if(out_press.gt.0.and.mod((itime-istart),oi_chan).eq.0)
     &          call output_press(out_press)
c
c     compute the divergence and the cfl condition 
           if (mod((itime-istart),oi_cfl).eq.0) then
              call divg
              call courant(dt0,u,v,w)
              write(26,*)' '
c
              if(FLSTLE.and.(.not.FLGVARY))then
c     recall ws is magnitude (>0)
c     recall gz has its own sign (gz<0)
c                 igmod=sqrt(gy(1,1)**2+gz(1,1)**2)
c                 igmod=1.d0/igmod
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k)
                 do k=1,nz
                    do j=1,ny
                       do i=1,nx
                          rv(i,j,k)=v(i,j,k)+ws
     &                             *gy(1,k)/sqrt(gy(1,k)**2+gz(1,k)**2)
                          rw(i,j,k)=w(i,j,k)+ws
     &                             *gz(1,k)/sqrt(gy(1,k)**2+gz(1,k)**2)
                       enddo
                    enddo
                 enddo                 
C$OMP END PARALLEL DO
                 write(26,*)'  Disperse phase info:'
                 write( *,*)'  Disperse phase info:'
                 call courant(dt0,u,rv,rw)
                 write(26,*)' '
              endif
           endif
c
c     output deposition and erosion
           if(.false..and.mod((itime-istart),oi_mean).eq.0)
     &       call output_2d_section('deposit',deposit,erosion,deposit,2)
c
c     output Chebyshev spectrum
           if (mod((itime-istart),oi_spec).eq.0)then
              call make_chb_spec(u,espz)
              call chb_spec_out(espz,'chspu',itime)
              call make_chb_spec(v,espz)
              call chb_spec_out(espz,'chspv',itime)
              call make_chb_spec(w,espz)
              call chb_spec_out(espz,'chspw',itime)
              if(FLTHRM)then
                 call make_chb_spec(tt,espz)
                 call chb_spec_out(espz,'chspT',itime)
              endif
           endif
c
c     output 1d stats
           if (.false..and.mod((itime-istart),oi_1d).eq.0) then
              call make_output_1d(30)
              n1d = n1d + 1
           endif
c     
c     output global balances
           if (.false..and.mod((itime-istart),oi_gbal).eq.0) then
              call make_gbalrhs(KE,Ep,diss,ddtEp)
              write(31,200) rtime,KE,Ep,diss,ddtEp
              ngbal = ngbal + 1
           endif
c
c     output mean profiles for some time steps
c$$$           if(mod((itime-istart),oi_cfl).eq.0) then
c$$$              open(11,file='test-data.dat',status='unknown')
c$$$              do k=1,nz
c$$$                 write(11,11)zets(k),v(nxh,nyh,k),vm(k,1)/dble(NT(1))
c$$$     &                ,vzm(k,1)/dble(NT(1))/re
c$$$     &                ,vwm(k,1)/dble(NT(1))
c$$$     &                -vm(k,1)/dble(NT(1))*wm(k,1)/dble(NT(1))
c$$$     &                ,tt(nxh,nyh,k),ttm(k,1)/dble(NT(1))
c$$$     &                ,ttzm(k,1)/dble(NT(1))/re/pr
c$$$     &                ,wttm(k,1)/dble(NT(1))
c$$$     &                -wm(k,1)/dble(NT(1))*ttm(k,1)/dble(NT(1))
c$$$              enddo
c$$$              close(11)
c$$$           end if
c 11   FORMAT(g25.15,g25.15,g25.15,g25.15,g25.15
c     &          ,g25.15,g25.15,g25.15,g25.15,g25.15)
 11   FORMAT(10(3x,g15.8))
c
c     transform back to Fourier space
           call fft_r2f_2d_new(u)
           call fft_r2f_2d_new(v)
           call fft_r2f_2d_new(w)
           if (FLTHRM) call fft_r2f_2d_new(tt)
c     
        endif!outputs in real space
c
c     close and reopen gbal and 1dout when main data written, 
c     so no data is lost during a crash and restart
        if (.false..and.mod((itime-istart),oi_chan).eq.0) then
           close(30)
           open(30,file='1dout')
           rewind(30)
           do j=1,n1d
c check in make_output_1d # of variables to read
              read(30,*) diss,diss,diss,diss,diss
           enddo
           close(31)
           open(31,file='gbal')
           rewind(31)
           do j=1,ngbal
              read(31,*) diss,diss,diss,diss,diss
           enddo
        endif
c
c     if the time step is divisible by oi_timer, print timers
        if (FLTIMER.and.mod((itime-istart),oi_timer).eq.0) then
           print *, 'Intermediate timings, time = ', itime
           call print_timers
        endif
c
c     all done with one complete time step
      enddo!itime
ccccccccccccccccccccccccccccccccccccccccccccc
c      write(40,110) (zets(k),u(nxh,nyh,k),k=1,nz)
!      write(40,111) (zets(k),vm( k,1)/dble(NT(1)),
!     &                       ttm(k,1)/dble(NT(1)),k=1,nz)
c      write(40,111) (zets(k)+rlenz/2.d0,v(1,1,k),dsin(rtime)
c     &-exp(-(zets(k)+rlenz/2.d0)*sqrt(re/2.d0))*dsin(rtime
c     &-(zets(k)+rlenz/2.d0)*sqrt(re/2.d0)),k=1,nz)
c      write(40,111) (zets(k)+rlenz/2.d0,v(1,1,k),dsin(rtime)
c     &-exp(-(-zets(k)+rlenz/2.d0)*sqrt(re/2.d0))*dsin(rtime
c     &-(-zets(k)+rlenz/2.d0)*sqrt(re/2.d0)),k=1,nz)
ccccccccccccccccccccccccccccccccccccccccccccc
c
c     fill in and output remaining data for the final time-step
      call advance(1)
c
      close(26)
      close(30)
      close(31)
cccccccccccccccccccccccccccccccccccccccccc GA ccccccccccccccccc
      close(unit=1110)
      close(unit=1111)
      close(unit=1112)
      close(unit=1120)
      close(unit=1121)
      close(unit=1122)
      close(unit=114)
      close(unit=115)
      close(unit=116)
      close(unit=120)
      close(unit=121)
      close(unit=122)
      close(unit=123)
      close(unit=124)
      close(unit=125)
cccccccccccccccccccccccccccccccccccccccccGA

c ================================================== c      
c      call finalize_ellip
c ================================================== c      
c
#ifdef FFTW3	
      call cleanfft
#endif
c
      if (FLTIMER) then
        print *, 'Final timings'
        call print_timers
      endif
c
 101  FORMAT(7i5,3f10.5)
 102  FORMAT(3f10.5,2i5)
 103  FORMAT(4i5) 
 109  FORMAT(i5,2f10.5)
 110  FORMAT(f12.8,5x,f12.8)
 111  FORMAT(f12.8,5x,f12.8,5x,f12.8)
 200  FORMAT(1p5e18.10)

      end subroutine



      PROGRAM main
      use common_m
      call common_allocate
      call main_run
      end
