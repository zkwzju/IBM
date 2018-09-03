module common_m
  use global_m
  use fft_m
  use flags_m
  use timers_m
  IMPLICIT none

  INTEGER istart,itime
  INTEGER oi_chan,out_press
  INTEGER oi_timer,oi_cfl,oi_spec,oi_mean,oi_gbal,oi_1d
  INTEGER msvx, msvy
  real*8, parameter :: pi = 3.14159265358979324d0
  REAL*8 re,rlenx,rleny,rlenz,hrlenz,deltax,deltay,idx,idy
  REAL*8 dt,dt_old,rtime,ipin,cflmax,gradpbar,ppA,ppT,cfpg,theta
  REAL*8 pr,ra,Tmax 
  REAL*8 uvab,   uvat,   uvbb,   uvbt,    ugb,    ugt,vgb,vgt
  REAL*8  wab,    wat,    wbb,    wbt,    wgb,    wgt
  REAL*8 ws,rey_p
  REAL*8 cferos1,cferos2 
  REAL*8 cfnl1v(3),cfnl2v(3),cfdifv(3),cfdelv(3),cfcumv(3)
  REAL*8 epssvx,epssvy
  REAL*8 uvmbct1,uvmbcb1,uvmbctn,uvmbcbn
  REAL*8 wmbct1,wmbcb1,wmbctn,wmbcbn
  REAL*8 ttmbct1,ttmbcb1,ttmbctn,ttmbcbn
  INTEGER kts(11)
  INTEGER, allocatable,dimension(:) :: NT
  integer irkk,n_l_max
  INTEGER flagibm,ibm_moving,icount
  INTEGER, allocatable,dimension(:) :: n_l


  CHARACTER*32 ch_file,tt_file,format_mode,sp_file
  real*8 volibm,volibm_mid
  real*8  ibmb(3),FI
  real*8,parameter,dimension(3) :: ibme=(/1.0d0,0.d0,0.d0/)
  real*8,allocatable,dimension(:) :: ttab,ttat,ttbb,ttbt,ttgb,ttgt
  real*8,allocatable,dimension(:) :: zets,deltaz,idz,ndeltaz,nidz,wzet
!!!!!!!!!!!!!!!!!!!!!

  real*8,allocatable,dimension(:) :: wavex,wavexs,wavey,waveys
  real*8,allocatable,dimension(:,:) :: chb
  real*8,allocatable,dimension(:,:) :: chbd1,chbd2
  real*8,allocatable,dimension(:,:) :: cn_drdr,ci_drdr
  real*8,allocatable,dimension(:) :: e_drdr
  real*8,allocatable,dimension(:,:) :: uvcn_mimi,uvci_mimi
  real*8,allocatable,dimension(:) :: uve_mimi
  real*8,allocatable,dimension(:,:) :: wcn_mimi,wci_mimi
  real*8,allocatable,dimension(:) :: we_mimi
  real*8,allocatable,dimension(:,:) :: ttcn_mimi,ttci_mimi
  real*8,allocatable,dimension(:) :: tte_mimi
  real*8,allocatable,dimension(:,:) :: cn_nodr,ci_nodr
  real*8,allocatable,dimension(:) :: e_nodr
  real*8,allocatable,dimension(:,:) :: cn_nono,ci_nono
  real*8,allocatable,dimension(:) :: e_nono
  real*8,allocatable,dimension(:,:,:) :: mlt
  real*8,allocatable,dimension(:,:) :: pbxnb,pbynb
  real*8,allocatable,dimension(:,:) :: pbxnt,pbynt
  real*8,allocatable,dimension(:,:) :: pbxob,pbyob
  real*8,allocatable,dimension(:,:) :: pbxot,pbyot
  real*8,allocatable,dimension(:,:) :: pcxnb,pcynb
  real*8,allocatable,dimension(:,:) :: pcxnt,pcynt
  real*8,allocatable,dimension(:,:) :: pcxob,pcyob
  real*8,allocatable,dimension(:,:) :: pcxot,pcyot

  real*8,allocatable,dimension(:) ::   prbc1,  prbcn,uvprbc1, uvprbcn,wprbc1, wprbcn,ttprbc1,ttprbcn
  real*8,allocatable,dimension(:,:,:) :: u,v,w,p,u_temp,v_temp,w_temp
  real*8,allocatable,dimension(:,:,:) :: ru,rv,rw,rp
  real*8,allocatable,dimension(:,:,:) :: h1,h2,h3
  real*8,allocatable,dimension(:,:,:) :: auxt1,auxt2
  real*8,allocatable,dimension(:,:,:) :: tt,rtt,htt
  real*8,allocatable,dimension(:,:,:) :: ul,vl,wl
  real*8,allocatable,dimension(:,:,:) :: ddxul,ddyul,ddzul
  real*8,allocatable,dimension(:,:) :: diff_u,diff_v,diff_w
  real*8,allocatable,dimension(:) :: svx,svy
  real*8,allocatable,dimension(:,:) :: deposit,erosion
  real*8,allocatable,dimension(:,:) :: ubct,ubcb,ubctax,ubcbax
  real*8,allocatable,dimension(:,:) :: vbct,vbcb,vbctax,vbcbax
  real*8,allocatable,dimension(:,:) :: wbct,wbcb
  real*8,allocatable,dimension(:,:) :: ttbct,ttbcb
  real*8,allocatable,dimension(:,:) :: gy,gz
  real*8,allocatable,dimension(:,:) ::  um,vm,wm,ttm&
       ,u2m,v2m,w2m,tt2m&
       ,u3m,v3m,w3m,tt3m&
       ,u4m,v4m,w4m,tt4m&
       ,uvm,uwm,uttm&
       ,vwm,vttm,wttm
  real*8,allocatable,dimension(:,:) ::  uxm,uym,uzm &
       ,ux2m,uy2m,uz2m&
       ,ux3m,uy3m,uz3m&
       ,ux4m,uy4m,uz4m
  real*8,allocatable,dimension(:,:) ::  vxm,vym,vzm&
       ,vx2m,vy2m,vz2m&
       ,vx3m,vy3m,vz3m&
       ,vx4m,vy4m,vz4m
  real*8,allocatable,dimension(:,:) ::  wxm,wym,wzm&
       ,wx2m,wy2m,wz2m&
       ,wx3m,wy3m,wz3m&
       ,wx4m,wy4m,wz4m
  real*8,allocatable,dimension(:,:) ::  ttxm,ttym,ttzm&
       ,ttx2m,tty2m,ttz2m&
       ,ttx3m,tty3m,ttz3m&
       ,ttx4m,tty4m,ttz4m
  real*8,allocatable,dimension(:,:,:) ::  uespxm,uespym  
  real*8,allocatable,dimension(:,:,:) ::  vespxm,vespym  
  real*8,allocatable,dimension(:,:,:) ::  wespxm,wespym  
  real*8,allocatable,dimension(:,:) ::  pwm,pm,u2wm,v2wm,uuzm,uwxm,vvzm,vwym,wwzm

  !     common variables used by the ibm
  real*8,allocatable,dimension(:) ::  xets,yets,xets4,yets4

  real*8,allocatable,dimension(:,:,:) ::  forcing_x,forcing_y,forcing_z,ddf_dum
  real*8,allocatable,dimension(:,:,:) ::  forcing_x1,forcing_y1,forcing_z1
  real*8,allocatable,dimension(:,:,:) ::  forcing_x2,forcing_y2,forcing_z2
  real*8,allocatable,dimension(:,:,:) ::  forcing_x3,forcing_y3,forcing_z3

  real*8,allocatable,dimension(:,:,:) :: forcing_x01,forcing_y01,forcing_z01
  real*8,allocatable,dimension(:,:,:) :: forcing_x02,forcing_y02,forcing_z02
  real*8,allocatable,dimension(:,:,:) :: forcing_x03,forcing_y03,forcing_z03
  
  real*8,allocatable,dimension(:,:,:,:,:) ::  ddf

  real*8,allocatable,dimension(:) ::  x_c,y_c,z_c,u_c,v_c,w_c,om_x,om_y,om_z
  real*8,allocatable,dimension(:) :: cell_v
  real*8,allocatable,dimension(:,:) ::  x_o,y_o,z_o,dv_l,ds
  real*8,allocatable,dimension(:,:) ::  rx_l,ry_l,rz_l
  real*8,allocatable,dimension(:) ::  c_d_total,c_lx_total,c_lz_total
  real*8,allocatable,dimension(:) ::  c_d_avg,c_lx_avg,c_lz_avg
  real*8,allocatable,dimension(:) ::  moment_x,moment_y,moment_z
  integer,allocatable,dimension(:,:) ::  p_iw,p_ie,p_js,p_jn,p_kb,p_kt
  real*8,allocatable,dimension(:,:) ::  u_p,v_p,w_p
  real*8,allocatable,dimension(:,:,:) ::  dpx,dpy,dpz
  real*8,allocatable,dimension(:) ::  pla_vol_fract
  real*8,allocatable,dimension(:,:,:) :: if_ibm

  real*8,allocatable,dimension(:) :: int_u0,int_v0,int_w0,&
       int_r_u0,int_r_v0,int_r_w0,&
       int_u1,int_v1,int_w1,&
       int_r_u1,int_r_v1,int_r_w1

  real*8,allocatable,dimension(:) :: int_fx,int_fy,int_fz
  real*8,allocatable,dimension(:) :: int_Fpx,int_Fpy,int_Fpz
  
contains
  subroutine common_allocate
    ! initial prameters
    call global_init
    ! fftw
    call fft_init

    !
    allocate( ttab(m),ttat(m),ttbb(m),ttbt(m),ttgb(m),ttgt(m) )
    allocate( zets(nz),deltaz(nzm),idz(nzm),ndeltaz(nzm),nidz(nzm),wzet(nz) )

    allocate( wavex(nx),wavexs(nx),wavey(ny),waveys(ny) )
    allocate( chb(nz0,nz) )
    allocate( chbd1(nz0,nz),chbd2(nz0,nz) )
    allocate( cn_drdr(nz0,nz),ci_drdr(nz0,nz),e_drdr(nz) )
    allocate( uvcn_mimi(nzmm0,nzmm),uvci_mimi(nzmm0,nzmm),uve_mimi(nzmm) )
    allocate( wcn_mimi(nzmm0,nzmm),wci_mimi(nzmm0,nzmm),we_mimi(nzmm) )
    allocate( ttcn_mimi(nzmm0,nzmm),ttci_mimi(nzmm0,nzmm),tte_mimi(nzmm) )
    allocate( cn_nodr(nz0,nz),ci_nodr(nz0,nz),e_nodr(nz) )
    allocate( cn_nono(nzmm0,nzmm),ci_nono(nzmm0,nzmm),e_nono(nzmm) )
    allocate( mlt(nz0,nx0,nyh) )
    allocate( pbxnb(nx0,nyh),pbynb(nx0,nyh) )
    allocate( pbxnt(nx0,nyh),pbynt(nx0,nyh) )
    allocate( pbxob(nx0,nyh),pbyob(nx0,nyh) )
    allocate( pbxot(nx0,nyh),pbyot(nx0,nyh) )
    allocate( pcxnb(nx0,nyh),pcynb(nx0,nyh) )
    allocate( pcxnt(nx0,nyh),pcynt(nx0,nyh) )
    allocate( pcxob(nx0,nyh),pcyob(nx0,nyh) )
    allocate( pcxot(nx0,nyh),pcyot(nx0,nyh) )

    allocate( prbc1(nzmm),  prbcn(nzmm),uvprbc1(nzmm), uvprbcn(nzmm),wprbc1(nzmm), wprbcn(nzmm),ttprbc1(nzmm),ttprbcn(nzmm) )
    allocate( u(nx0,ny0,nz),v(nx0,ny0,nz),w(nx0,ny0,nz),p(nx0,ny0,nz) )
    allocate( u_temp(nx0,ny0,nz),v_temp(nx0,ny0,nz),w_temp(nx0,ny0,nz) )
    allocate( ru(nx0,ny0,nz),rv(nx0,ny0,nz),rw(nx0,ny0,nz),rp(nx0,ny0,nz) )
    allocate( h1(nx0,ny0,nz),h2(nx0,ny0,nz),h3(nx0,ny0,nz) )

    allocate( auxt1(nz0,nx0,nyh),auxt2(nz0,nx0,nyh) )
    allocate( tt(nx0,ny0,nz),rtt(nx0,ny0,nz),htt(nx0,ny0,nz) )
    allocate( ul(nxl0,nyl0,nz),vl(nxl0,nyl0,nz),wl(nxl0,nyl0,nz) )
    allocate( ddxul(nxl0,nyl0,nz),ddyul(nxl0,nyl0,nz),ddzul(nxl0,nyl0,nz) )
    allocate( svx(nx),svy(nyh) )
    allocate( deposit(nx0,ny0),erosion(nx0,ny0) )
    allocate( ubct(nx0,ny0),ubcb(nx0,ny0),ubctax(nx0,ny0),ubcbax(nx0,ny0) )
    allocate( vbct(nx0,ny0),vbcb(nx0,ny0),vbctax(nx0,ny0),vbcbax(nx0,ny0) )
    allocate( wbct(nx0,ny0),wbcb(nx0,ny0) )
    allocate( ttbct(nx0,ny0),ttbcb(nx0,ny0) )
    allocate( gy(nyl0,nz),gz(nyl0,nz) )



    allocate(   um(nz,np),vm(nz,np),wm(nz,np),ttm(nz,np)&
         ,u2m(nz,np),v2m(nz,np),w2m(nz,np),tt2m(nz,np)&
         ,u3m(nz,np),v3m(nz,np),w3m(nz,np),tt3m(nz,np)&
         ,u4m(nz,np),v4m(nz,np),w4m(nz,np),tt4m(nz,np)&
         ,uvm(nz,np),uwm(nz,np),uttm(nz,np)&
         ,vwm(nz,np),vttm(nz,np)&
         ,wttm(nz,np) )
    allocate(   uxm( nz,np),uym( nz,np),uzm( nz,np) &
         ,ux2m(nz,np),uy2m(nz,np),uz2m(nz,np)&
         ,ux3m(nz,np),uy3m(nz,np),uz3m(nz,np)&
         ,ux4m(nz,np),uy4m(nz,np),uz4m(nz,np) )
    allocate(   vxm( nz,np),vym( nz,np),vzm( nz,np)&
         ,vx2m(nz,np),vy2m(nz,np),vz2m(nz,np)&
         ,vx3m(nz,np),vy3m(nz,np),vz3m(nz,np)&
         ,vx4m(nz,np),vy4m(nz,np),vz4m(nz,np) )
    allocate(   wxm( nz,np),wym( nz,np),wzm( nz,np)&
         ,wx2m(nz,np),wy2m(nz,np),wz2m(nz,np)&
         ,wx3m(nz,np),wy3m(nz,np),wz3m(nz,np)&
         ,wx4m(nz,np),wy4m(nz,np),wz4m(nz,np) )
    allocate(   ttxm( nz,np),ttym( nz,np),ttzm( nz,np)&
         ,ttx2m(nz,np),tty2m(nz,np),ttz2m(nz,np)&
         ,ttx3m(nz,np),tty3m(nz,np),ttz3m(nz,np)&
         ,ttx4m(nz,np),tty4m(nz,np),ttz4m(nz,np) )
    allocate(   uespxm(nxh,nz,np),uespym(nyh,nz,np) )
    allocate(   vespxm(nxh,nz,np),vespym(nyh,nz,np) )  
    allocate(   wespxm(nxh,nz,np),wespym(nyh,nz,np) )  
    allocate(   pwm(nz,np),pm(nz,np)&
         ,u2wm(nz,np),v2wm(nz,np)&
         ,uuzm(nz,np),uwxm(nz,np)&
         ,vvzm(nz,np),vwym(nz,np)&
         ,wwzm(nz,np) )

    !     common variables used by the ibm
    allocate(   xets(nx0),yets(ny0),xets4(-1:nx0+1),yets4(-1:ny0+1) )
    allocate(   forcing_x(nx0,ny0,nz),forcing_y(nx0,ny0,nz),forcing_z(nx0,ny0,nz) )
    
    allocate(   forcing_x1(nx0,ny0,nz),forcing_y1(nx0,ny0,nz),forcing_z1(nx0,ny0,nz) )
    allocate(   forcing_x2(nx0,ny0,nz),forcing_y2(nx0,ny0,nz),forcing_z2(nx0,ny0,nz) )
    allocate(   forcing_x3(nx0,ny0,nz),forcing_y3(nx0,ny0,nz),forcing_z3(nx0,ny0,nz) )
    allocate(   forcing_x01(nx0,ny0,nz),forcing_y01(nx0,ny0,nz),forcing_z01(nx0,ny0,nz) )
    allocate(   forcing_x02(nx0,ny0,nz),forcing_y02(nx0,ny0,nz),forcing_z02(nx0,ny0,nz) )
    allocate(   forcing_x03(nx0,ny0,nz),forcing_y03(nx0,ny0,nz),forcing_z03(nx0,ny0,nz) )
    
    allocate(   ddf(n_ll,num_p,3,3,4),ddf_dum(nx,ny,nz) )

    allocate(   x_c(num_p),y_c(num_p),z_c(num_p)&
         ,u_c(num_p),v_c(num_p),w_c(num_p)&
         ,om_x(num_p),om_y(num_p),om_z(num_p) )
    allocate( cell_v(nz) )    ! Eulerian cell volume
    allocate(   x_o(n_ll,num_p),y_o(n_ll,num_p),z_o(n_ll,num_p),dv_l(n_ll,num_p),ds(n_ll,num_p) )
    allocate(   rx_l(n_ll,num_p),ry_l(n_ll,num_p),rz_l(n_ll,num_p) )
    allocate(   c_d_total(num_p),c_lx_total(num_p),c_lz_total(num_p) )
    allocate(   c_d_avg(num_p),c_lx_avg(num_p),c_lz_avg(num_p) )
    allocate(   moment_x(num_p),moment_y(num_p),moment_z(num_p) )
    allocate(   p_iw(n_ll,num_p),p_ie(n_ll,num_p)&
         ,p_js(n_ll,num_p),p_jn(n_ll,num_p)&
         ,p_kb(n_ll,num_p),p_kt(n_ll,num_p) )
    allocate(   u_p(n_ll,num_p),v_p(n_ll,num_p),w_p(n_ll,num_p) )
    allocate(   dpx(nx0,ny0,nz),dpy(nx0,ny0,nz),dpz(nx0,ny0,nz) )
    allocate(   pla_vol_fract(nz),if_ibm(nx,ny,nz) )

    allocate( n_l(num_p), NT(np) )

    allocate( int_u0(num_p),int_v0(num_p),int_w0(num_p),&
       int_r_u0(num_p),int_r_v0(num_p),int_r_w0(num_p),&
       int_u1(num_p),int_v1(num_p),int_w1(num_p),&
       int_r_u1(num_p),int_r_v1(num_p),int_r_w1(num_p) )

    allocate( int_fx(num_p),int_fy(num_p),int_fz(num_p))
    allocate( int_Fpx(num_p),int_Fpy(num_p),int_Fpz(num_p))

    allocate( diff_u(n_ll,num_p),diff_v(n_ll,num_p),diff_w(n_ll,num_p))
  end subroutine common_allocate


end module common_m
