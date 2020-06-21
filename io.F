      SUBROUTINE input
      use common_m

      implicit none
c
      INTEGER itold
      REAL*8  rtimet,dt_oldt
      INTEGER idum
      REAL*8  dum
      CHARACTER*32 anum

      open(11,file=anum(ch_file,istart),form=format_mode)
      rewind(11)
      read(11) rtime,itold,dt_old,idum,idum,idum,dum,dum,dum
      read(11) u,v,w
c      read(11) pbxnt,pcxnt,pbynt,pcynt,pbxot,pcxot,pbyot,pcyot
c      read(11) pbxnb,pcxnb,pbynb,pcynb,pbxob,pcxob,pbyob,pcyob
      close(11)
      write(26,*)'input: rtime,itold,dt_old ', rtime,itold,dt_old
      if (itold.ne.istart) then
        write(26,*)'Error in input: itold is not istart'
        stop
      endif
c
      if (FLTHRM) then
        open(11,file=anum(tt_file,istart),form=format_mode)
        rewind(11)
        read(11) rtimet,itold,dt_oldt,idum,idum,idum,dum,dum,dum
        read(11) tt
        close(11)
        if (itold.ne.istart) then
          write(26,*)'Error in input: itold is not istart'
          stop
        endif
        if ((rtime.ne.rtimet).or.(dt_old.ne.dt_oldt)) then
          write(26,*)'Error in input: time mismatch with thermal file',
     &      rtime, rtimet, dt_old, dt_oldt
          stop
        endif
      endif
c
c if ws.ne.0 input deposit
      if(.false..and.FLSTLE)then 
         open(11,file=anum('deposit',istart),form=format_mode)
         rewind(11)
         read(11) rtimet,itold,dt_oldt,idum,idum,idum,dum,dum,dum,idum
         read(11) deposit,erosion
         close(11)
      endif
c
      return
      END
c
      SUBROUTINE output
      use common_m
      use ellip_common_m, only: lcheck,lunformatted,for_px,for_py,
     &     for_pz,torq_x,torq_y,torq_z,xlp,ylp,zlp,I_particle,
     &     rho_particle,V_particle,axis_1,axis_2,axis_3
      implicit none

c
      CHARACTER*32 anum
      INTEGER      n,j
c
      if (FLTIMER) call begin_timer(30)
      open(21,file=anum(ch_file,itime),form=format_mode)
      write(21) rtime,itime,dt,nx0,ny0,nz,deltax,deltay,rlenz
      write(21) u,v,w
      write(21) pbxnt,pcxnt,pbynt,pcynt,pbxot,pcxot,pbyot,pcyot
      write(21) pbxnb,pcxnb,pbynb,pcynb,pbxob,pcxob,pbyob,pcyob
      close(21)
c     
      if(FLTHRM)then
         open(21,file=anum(tt_file,itime),form=format_mode)
         write(21) rtime,itime,dt,nx0,ny0,nz,deltax,deltay,rlenz
         write(21) tt
         close(21)
      endif
cccccccccccccccccccccccccc GA cccccccccccccccccccccccccccccccccc
      open(126,file=anum(sp_file,itime),form=format_mode)
      write(126) rho_particle,V_particle,axis_1,axis_2,axis_3,
     &           I_particle,x_c,y_c,z_c,u_c,v_c,
     &           w_c,om_x,om_y,om_z,
     &           for_px,for_py,for_pz,
     &           torq_x,torq_y,torq_z,
     &           n_l,xlp,ylp,zlp
      close(126)
cccccccccccccccccccccccccc GA
c

      if(lcheck) then
         open(375,file=anum('fIBM_1',itime),form='unformatted')
         write(375) forcing_x1,forcing_y1,forcing_z1
         close(375)
         
         open(375,file=anum('fIBM_2',itime),form='unformatted')
         write(375) forcing_x2,forcing_y2,forcing_z2
         close(375)
         
         open(375,file=anum('fIBM_3',itime),form='unformatted')
         write(375) forcing_x3,forcing_y3,forcing_z3
         close(375)

         open(376,file=anum('fIBM_01',itime),form='unformatted')
         write(376) forcing_x01,forcing_y01,forcing_z01
         close(376)

         open(376,file=anum('fIBM_02',itime),form='unformatted')
         write(376) forcing_x02,forcing_y02,forcing_z02
         close(376)
         
         open(376,file=anum('fIBM_03',itime),form='unformatted')
         write(376) forcing_x03,forcing_y03,forcing_z03
         close(376)
 
      endif

      
      if (FLTIMER) call end_timer(30)
      return
      END
cccccccccccccccccccccccccccccccccccccccccccc GA ccccccccccccccccccccc
      SUBROUTINE output_sphere
      use common_m

c
      CHARACTER*32 anum
c
      open(21,file=anum(ch_file,itime),form=format_mode)
      write(21) rtime,itime,dt,nx0,ny0,nz,deltax,deltay,rlenz
      write(21) u,v,w
      write(21) pbxnt,pcxnt,pbynt,pcynt,pbxot,pcxot,pbyot,pcyot
      write(21) pbxnb,pcxnb,pbynb,pcynb,pbxob,pcxob,pbyob,pcyob
      close(21)
c
      if(FLTHRM)then
         open(21,file=anum(tt_file,itime),form=format_mode)
         write(21) rtime,itime,dt,nx0,ny0,nz,deltax,deltay,rlenz
         write(21) tt
         close(21)
      endif
c
      return
      END
cccccccccccccccccccccccccccccccccccccccccccc GA 
c
      SUBROUTINE output_2d_section(name,val2d_1,val2d_2,val2d_3,dim)
      use common_m

      REAL*8 val2d_1(nx0,ny0),val2d_2(nx0,ny0),val2d_3(nx0,ny0)
      INTEGER dim
      CHARACTER*(*) name
c
      CHARACTER*32 anum
c
      if (FLTIMER) call begin_timer(36)
c write averaged value
      open(21,file=anum(name,itime),form=format_mode)
      write(21) rtime,itime,dt,nx0,ny0,nz,deltax,deltay,rlenz,0
      if(dim.eq.1) write(21) val2d_1
      if(dim.eq.2) write(21) val2d_1,val2d_2
      if(dim.eq.3) write(21) val2d_1,val2d_2,val2d_3
      close(21)
c
      if (FLTIMER) call end_timer(36)
      return
      END
c
      SUBROUTINE input_means
      use common_m

c
      CHARACTER*32 anum
      INTEGER nz_old,np_old,idum
      REAL*8 rlenzaux,dum
c
c      if (FLTIMER) call begin_timer()
      open(21,file=anum('means',itime),form=format_mode)
      read(21)dum,idum,dum,idum,idum,nz_old,np_old,dum,dum,rlenzaux,NT
      read(21)dum,dum,dum
      if(nz.ne.nz_old.or.np.ne.np_old.or.rlenz.ne.rlenzaux)then
         write(26,*)'Error in input of mean: nz_old ne nz'
         write(26,*)'                        or rlenz_old ne rlenz'
         write(26,*)'                        or np_old ne np'
         stop
      endif
      read(21)um,vm,wm,ttm,u2m,v2m,w2m,tt2m,u3m,v3m,w3m,tt3m,u4m,v4m
     &,w4m,tt4m,uvm,uwm,uttm,vwm,vttm,wttm,uxm,uym,uzm,ux2m,uy2m,uz2m
     &,ux3m,uy3m,uz3m,ux4m,uy4m,uz4m,vxm,vym,vzm,vx2m,vy2m,vz2m,vx3m
     &,vy3m,vz3m,vx4m,vy4m,vz4m,wxm,wym,wzm,wx2m,wy2m,wz2m,wx3m,wy3m
     &,wz3m,wx4m,wy4m,wz4m,ttxm,ttym,ttzm,ttx2m,tty2m,ttz2m,ttx3m,tty3m
     &,ttz3m,ttx4m,tty4m,ttz4m
      read(21)uespxm,uespym,vespxm,vespym,wespxm,wespym
      read(21)pwm,pm,u2wm,v2wm,uuzm,uwxm,vvzm,vwym,wwzm
      close(21)
c
c      if (FLTIMER) call end_timer()
      return
      END
c
      SUBROUTINE output_means
      use common_m

c
      CHARACTER*32 anum
c
c      if (FLTIMER) call begin_timer()
      open(21,file=anum('means',itime),form=format_mode)
      write(21)rtime,itime,dt,nx0,ny0,nz,np,deltax,deltay,rlenz,NT
      write(21)re,pr,ra
      write(21)um,vm,wm,ttm,u2m,v2m,w2m,tt2m,u3m,v3m,w3m,tt3m,u4m,v4m
     &,w4m,tt4m,uvm,uwm,uttm,vwm,vttm,wttm,uxm,uym,uzm,ux2m,uy2m,uz2m
     &,ux3m,uy3m,uz3m,ux4m,uy4m,uz4m,vxm,vym,vzm,vx2m,vy2m,vz2m,vx3m
     &,vy3m,vz3m,vx4m,vy4m,vz4m,wxm,wym,wzm,wx2m,wy2m,wz2m,wx3m,wy3m
     &,wz3m,wx4m,wy4m,wz4m,ttxm,ttym,ttzm,ttx2m,tty2m,ttz2m,ttx3m,tty3m
     &,ttz3m,ttx4m,tty4m,ttz4m
      write(21)uespxm,uespym,vespxm,vespym,wespxm,wespym
      write(21)pwm,pm,u2wm,v2wm,uuzm,uwxm,vvzm,vwym,wwzm
      close(21)
c
c      if (FLTIMER) call end_timer()
      return
      END
c
      SUBROUTINE mean5_out(cl1,cl2,cl3,cl4,cl5,m1,m2,m3,m4,m5,outfile,t)
      use common_m

      INTEGER t
      REAL*8 cl1(nz),cl2(nz),cl3(nz),cl4(nz),cl5(nz),m1,m2,m3,m4,m5
      CHARACTER*(*) outfile
c
      INTEGER k
      CHARACTER*32 anum
c
      open(21,file=anum(outfile,t))
      write(21,1000) m1, m2, m3, m4, m5
      do k=1,nz
        write(21,2000) zets(k), cl1(k), cl2(k), cl3(k), cl4(k), cl5(k)
      enddo
      close(21)
c
      return
 1000 FORMAT('#   means:        ', 1p9e18.10)
 2000 FORMAT(1p9e18.10)
      END
c
      SUBROUTINE spec_out(espx,dspx,espy,dspy,outfile,t)
      use common_m

      INTEGER t
      REAL*8 espx(nxh,nz),dspx(nxh,nz),espy(nyh,nz),dspy(nyh,nz)
      CHARACTER*(*) outfile
c
      CHARACTER*32 anum
c
      open(21,file=anum(outfile,t),form=format_mode)
      write(21) rtime,itime,dt,nx0,ny0,nz,deltax,deltay,rlenz
      write(21) espx
      write(21) dspx
      write(21) espy
      write(21) dspy
      close(21)
c
      return
      END
c
      SUBROUTINE chb_spec_out(espz,outfile,t)
      use common_m

      INTEGER t
      REAL*8 espz(nz)
      CHARACTER*(*) outfile
c
      CHARACTER*32 anum
c
      open(21,file=anum(outfile,t),form=format_mode)
      write(21) rtime,itime,dt,nx0,ny0,nz,deltax,deltay,rlenz
      write(21) espz
      close(21)
c
      return
      END
c
      SUBROUTINE output_press(code)
      use common_m

      INTEGER code
c
c     This subroutine postprocesses and outputs pressure for output 
c     according to code:
c     code=1 -> output dyn press + hydro press due to strat by variable tt
c     code=2 -> output dyn press
c     Needs tt in physical space
c
      real*8 tmeanp,hydp,s
      integer i,j,k,l
      CHARACTER*32 anum
c
c      if (FLTIMER) call begin_timer(36)
c
c     copy p to rw
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(k)
      do k=1,nz
        call dcopy(nx0y,p(1,1,k),1,rw(1,1,k),1)
      enddo
C$OMP END PARALLEL DO
c
c     IFFT rw -> rw is pressure
      call fft_f2r_2d_new(rw)
c
c     substract stratification pressure 
c     and compute total mean dynamic pressure
      tmeanp = 0.d0
      if(code.eq.2)then
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k,l,hydp,s)
C$OMP&                           ,REDUCTION(+:tmeanp)
         do k=1,nz
            s=0.d0
            do j=1,ny
               do i=1,nx
                  hydp=0.d0
                  do l=1,k
                     hydp=hydp+tt(i,j,l)*(rlenz*wzet(l))
                  enddo
                  rw(i,j,k)=rw(i,j,k)-hydp
                  s=s+rw(i,j,k)
               enddo
            enddo
            s=s/dble(nx*ny)
            tmeanp=tmeanp+s*wzet(k)
         enddo
C$OMP END PARALLEL DO
      elseif(code.eq.1)then
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k,s)
C$OMP&                           ,REDUCTION(+:tmeanp)
         do k=1,nz
            s=0.d0
            do j=1,ny
               do i=1,nx
                  s=s+rw(i,j,k)
               enddo
            enddo
            s=s/dble(nx*ny)
            tmeanp=tmeanp+s*wzet(k)
         enddo
C$OMP END PARALLEL DO
      endif
c
c     substract mean pressure to normalize output
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k)
         do k=1,nz
            do j=1,ny
               do i=1,nx
                  rw(i,j,k)=rw(i,j,k)-tmeanp
               enddo
            enddo
         enddo
C$OMP END PARALLEL DO
c
      open(21,file=anum('press',itime),form=format_mode)
      write(21) rtime,itime,dt,nx0,ny0,nz,deltax,deltay,rlenz
      write(21) rw
      close(21)
c
c      if (FLTIMER) call end_timer(36)
      return
      END
c
      SUBROUTINE output_adcpdata
C     NOTE: this subroutine assumes physical space data in ru,rv,rw and rtt
      use common_m

      integer nxj,nyj
      parameter(nxj=21,nyj=21)
      real   ucil(nxj,nyj,nz),vcil(nxj,nyj,nz),wcil(nxj,nyj,nz)
      real   uslicel(3,ny,nz),vslicel(3,ny,nz),wslicel(3,ny,nz)
      real   uslicet(nx,3,nz),vslicet(nx,3,nz),wslicet(nx,3,nz)
      integer i,j,k
      CHARACTER*32 anum
ccccccccccccccccccccccccc
c      real   ucm(nz),vcm(nz),wcm(nz)
ccccccccccccccccccccccccc
c
      do k=1,nz
         do j=1,nyj
            do i=1,nxj
               ucil(i,j,k)=ru(nxhp-nxj/2-1+i,nyhp-nyj/2-1+j,k)
               vcil(i,j,k)=rv(nxhp-nxj/2-1+i,nyhp-nyj/2-1+j,k)
               wcil(i,j,k)=rw(nxhp-nxj/2-1+i,nyhp-nyj/2-1+j,k)
            enddo
         enddo
      enddo
ccccccccccccccccccccccccc
c      do k=1,nz
c         ucm(k)=0.0
c         vcm(k)=0.0
c         wcm(k)=0.0
c      enddo
c      do k=1,nz
c         do j=1,nyj
c            do i=1,nxj
c               ucm(k)=ucm(k)+ucil(i,j,k)/dble(nxj*nyj)
c               vcm(k)=vcm(k)+vcil(i,j,k)/dble(nxj*nyj)
c               wcm(k)=wcm(k)+wcil(i,j,k)/dble(nxj*nyj)
c            enddo
c         enddo
c      enddo
c      write(12,*)(zets(k),ucm(k),vcm(k),wcm(k),z=1,nz)
c      write(12,*)(zets(k),ucil(1,1,k),vcil(1,1,k),wcil(1,1,k),k=1,nz)
ccccccccccccccccccccccccc
c
      do k=1,nz
         do j=1,ny
            do i=1,3
               uslicel(i,j,k)=ru(nxhp-2+i,j,k)
               vslicel(i,j,k)=rv(nxhp-2+i,j,k)
               wslicel(i,j,k)=rw(nxhp-2+i,j,k)
            enddo
         enddo
      enddo
c
      do k=1,nz
         do j=1,3
            do i=1,nx
               uslicet(i,j,k)=ru(i,nyhp-2+j,k)
               vslicet(i,j,k)=rv(i,nyhp-2+j,k)
               wslicet(i,j,k)=rw(i,nyhp-2+j,k)
            enddo
         enddo
      enddo
c
      open(21,file=anum('vel-cil',itime),form=format_mode)
      write(21) nxj,nyj,nz,re,rlenx,rleny,rlenz
      write(21) ucil,vcil,wcil
      close(21)
c
      open(21,file=anum('vel-slicel',itime),form=format_mode)
      write(21) 3,ny,nz,re,rlenx,rleny,rlenz
      write(21) uslicel,vslicel,wslicel
      close(21)
c      
      open(21,file=anum('vel-slicet',itime),form=format_mode)
      write(21) nx,3,nz,re,rlenx,rleny,rlenz
      write(21) uslicet,vslicet,wslicet
      close(21)
c      
      return
      END      
c
      FUNCTION anum(string,num)
      IMPLICIT none
      CHARACTER*32 anum
      CHARACTER*(*) string
      INTEGER num
c
      INTEGER i,n,dot,fin
c
      anum = string
      dot = index(anum,' ')
      if (dot.gt.0.and.num.ge.0) then
        if (dot.le.32) anum(dot:dot) = '.'
        fin = dot
        n = num
 100    fin = fin + 1
        n = n/10
        if (n.gt.0) go to 100
        n = num
        do i=fin,dot+1,-1
          if (i.le.32) anum(i:i) = char(ichar('0')+mod(n,10))
          n = n/10
        enddo
      endif
c
      return
      END
