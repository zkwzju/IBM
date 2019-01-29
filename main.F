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

           endif

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


c     transform back to Fourier space
           call fft_r2f_2d_new(u)
           call fft_r2f_2d_new(v)
           call fft_r2f_2d_new(w)
           if (FLTHRM) call fft_r2f_2d_new(tt)
c     
        endif!outputs in real space
c

c     if the time step is divisible by oi_timer, print timers
        if (FLTIMER.and.mod((itime-istart),oi_timer).eq.0) then
           print *, 'Intermediate timings, time = ', itime
           call print_timers
        endif
c
c     all done with one complete time step
      enddo!itime

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
      close(unit=223)
      close(unit=224)
      close(unit=225)
      
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
