module ellip_common_m
  use precision_m
  use rng_m
  use common_m, only: x_o,y_o,z_o,x_c,y_c,z_c,rx_l,ry_l,rz_l,u_p,v_p,w_p,u_c,v_c,w_c,om_x,om_y,om_z,xets,yets,zets,p_iw,p_ie,p_js,p_jn,p_kb,p_kt,deltaz,deltax,deltay,rlenx,rleny,rlenz,hrlenz,n_l,n_ll,nx,ny,nz,num_p,dv_l,itime,u,v,w,cell_v,pi
  implicit none
  integer, allocatable :: x_0(:),x_1(:),y_0(:),y_1(:),z_0(:),z_1(:)
  real(wp), allocatable :: xlp(:),ylp(:),zlp(:)
  real(wp)              :: S_ellp   ! aproximate surface area
  real(wp)              :: vol_ellip
  real(wp)              :: om_ellip(3),om_ellip_b(3)
  real(wp)              :: I_ellip(3)
  real(wp), allocatable :: axis_1(:,:),axis_2(:,:),axis_3(:,:)
  real(wp)              :: rho_p,rho_f,vt_l,vb_l
  real(wp)              :: grav(3)
  real(wp), allocatable :: for_px(:),for_py(:),for_pz(:)
  real(wp), allocatable :: torq_x(:),torq_y(:),torq_z(:)

  real(wp) :: xa,xb,xc
  real(wp) :: mesh_length
  real(wp) :: v0_l
  real(wp) :: scale_dv,scale_p
  real(wp) :: frac_step,aalpha
  integer  :: imethod
  integer  :: iweight,iflag_ibm
  logical  :: lexact_weight,lupdate_dv,lrigid,lnonuniform,lrotation,ltranslation
  logical  :: lclip_x,lclip_y,lclip_z,lclip_ox,lclip_oy,lclip_oz
  logical  :: lcheck,llinear_v,lsphere,lconvect
  character(len=100) :: sname,sname_w
end module ellip_common_m
