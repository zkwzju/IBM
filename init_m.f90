module init_m
  use common_m
  INTEGER m_fin,sph_fin,iters,itfirst,iend,i,j,k
  INTEGER tt_yes,ngbal,n1d
  REAL*8  dt0,dt_p,KE,Ep,diss,ddtEp

contains
  subroutine initialization
    use parser_m
    use ellipsoid_m
    implicit none
    integer  :: reclen
    logical  :: lviscosity
    real(wp) :: A_size(3,3)
    character(len=120) :: file_position,file_form,file_status
    
    print *, 'Program running'


    
    !  ch_fin = 0 generates an initial solution
    !  ch_fin = 1 reads the initial data from ch_file.(istart)
    !  outfirst = whether to output data (other than c and u) at first time-step

    call parser_read('ch_file',ch_file,'vel')
    call parser_read('sp_file',sp_file,'sph')
    call parser_read('rlenx',rlenx,10.0_wp)
    call parser_read('rleny',rleny,10.0_wp)
    call parser_read('rlenz',rlenz,10.0_wp)
    call parser_read('uvbb',uvbb,0.0_wp)
    call parser_read('uvbt',uvbt,0.0_wp)
    call parser_read('istart',istart,0)
    call parser_read('iters',iters,1000)
    call parser_read('ch_fin',ch_fin,0)
    call parser_read('m_fin',m_fin,0)
    call parser_read('oi_chan',oi_chan,100.0_wp)
    call parser_read('out_press',out_press,100)
    call parser_read('oi_timer',oi_timer,1000000)
    call parser_read('oi_cfl',oi_cfl,100)
    call parser_read('oi_mean',oi_mean,1000000)
    call parser_read('oi_spec',oi_spec,1000000)
    call parser_read('oi_gbal',oi_gbal,1000000)
    call parser_read('oi_1d',oi_1d,1000000)
    call parser_read('Reynolds number',re,1.0_wp)
    call parser_is_defined('Kinematic viscosity',lviscosity)
    if(lviscosity) then
       call parser_read('Kinematic viscosity',re)
       re = 1.0_wp/re                         ! re is the inverse of viscosity
    endif
    call parser_read('Pressure gradient',ppA,0.0_wp)
    call parser_read('Time step dt',dt_p,0.1_wp)
    call parser_read('Max CFL',cflmax,1.0_wp)
    call parser_read('ppT',ppT,0.0_wp)
    call parser_read('Channel tilt angle',theta,0.d0)
    call parser_read('Including temperature',tt_yes,0)

    hrlenz = rlenz/2.0d0
    kts = (/5,10,30,60,90,120,130,140,150,160,180/)  

    dt = dt_p
    FLTIMER = oi_timer.gt.0
    FLAVER  = oi_mean.gt.0
    FLTHRM  = tt_yes.ne.0
    FLPGVAR = ppT.gt.0.d0

    if (FLTIMER) then
       call init_timers
    else
       oi_timer = 1
    endif

!!$    if (FLTHRM) then
!!$       read '(a)', tt_file
!!$       read(*,*) pr,ra,Tmax
!!$       read(*,*) ttab(1),ttat(1) 
!!$       read(*,*) ttbb(1),ttbt(1)
!!$       read(*,*) ttgb(1),ttgt(1)
!!$       read(*,*) ws,rey_p
!!$       read '(a)', identifier
!!$       if (identifier(1:9).ne.'#END_THRM') then
!!$          print *, 'Expecting identifier #END_THRM in input file'
!!$          stop
!!$       endif
!!$    endif

    FLSTLE  = .false.

    epssvx=0.d0
    epssvy=0.d0     
!!$    if (FLTHRM) then
!!$       read '(a)', identifier
!!$       if (identifier(1:8).ne.'#SP_VISC') then
!!$          print *, 'Expecting identifier #SP_VISC in input file'
!!$          stop
!!$       endif
!!$       read(*,*) epssvx,epssvy
!!$       read(*,*) msvx,msvy 
!!$       read '(a)', identifier
!!$       if (identifier(1:12).ne.'#END_SP_VISC') then
!!$          print *, 'Expecting identifier #END_SP_VISC in input file'
!!$          stop
!!$       endif
!!$    endif
    FLSV_NO=.false.
    FLSV_YES=epssvx.gt.0.or.epssvy.gt.0

    itfirst = istart+1
    iend = istart + iters
    itime = istart


    open(26,file='logfile')
    call initial(dt0,m_fin)  ! initialization of IBM related
    call correlation(n_ll)          ! to check weight if too big
    
    ! output grid
    open(100,file='grid.dat')
    write(100,*) 'x',nx+1
    write(100,*) xets(1:nx+1)
    write(100,*) 'y',ny+1
    write(100,*) yets(1:ny+1)
    write(100,*) 'z',nz
    write(100,*) zets(1:nz)
    close(100)


    
    ! produce complete output for initial condition if ...
    if(ch_fin.eq.0)then
       call output
       file_position = 'rewind'
       file_status   = 'unknown'
    else
       file_position = 'append'
       file_status   = 'old'
    endif

    call divg

    ! transform variables to Fourier space
    call fft_r2f_2d_new(u) 
    call fft_r2f_2d_new(v) 
    call fft_r2f_2d_new(w) 
    if (FLTHRM) call fft_r2f_2d_new(tt)

    ! transform bc's (that were computed in initial)
    ! NOTE: to do time varying bc see subroutine advance
    call fft_r2f_2d_new_slice(ubct)
    call fft_r2f_2d_new_slice(ubcb)
    call fft_r2f_2d_new_slice(vbct)
    call fft_r2f_2d_new_slice(vbcb)
    call fft_r2f_2d_new_slice(wbct)
    call fft_r2f_2d_new_slice(wbcb)
    if (FLTHRM) call fft_r2f_2d_new_slice(ttbct)
    if (FLTHRM) call fft_r2f_2d_new_slice(ttbcb)

    if(FLPGVAR)then
       cfpg=1.d0/ppT
    else
       gradpbar=ppA
    endif
    
    open(unit=114, file='liftx',status=file_status,position=file_position)
    open(unit=115, file='drag', status=file_status,position=file_position)
    open(unit=116, file='liftz',status=file_status,position=file_position)


    if(lunformatted) then
       file_form = 'unformatted'
    else
       file_form = 'formatted'
    endif
    
    if(lcheck) then
          ! open file for velocity difference between fluid and Lagrangian marker
       open(666,file='vel_diff.dat',status=file_status,position=file_position,form=file_form)
    endif
    open(unit=88,file='orientation_matrix',status=file_status,position=file_position,form=file_form)
    open(unit=1110, file='position_x'     ,status=file_status,position=file_position,form=file_form)
    open(unit=1111, file='position_y'     ,status=file_status,position=file_position,form=file_form)
    open(unit=1112, file='position_z'     ,status=file_status,position=file_position,form=file_form)
    open(unit=1120, file='velocity_x'     ,status=file_status,position=file_position,form=file_form)
    open(unit=1121, file='velocity_y'     ,status=file_status,position=file_position,form=file_form)
    open(unit=1122, file='velocity_z'     ,status=file_status,position=file_position,form=file_form)
    open(unit=120, file='om_x'            ,status=file_status,position=file_position,form=file_form)
    open(unit=121, file='om_y'            ,status=file_status,position=file_position,form=file_form)
    open(unit=122, file='om_z'            ,status=file_status,position=file_position,form=file_form)

    open(unit=123, file='moment_x'        ,status=file_status,position=file_position,form=file_form)
    open(unit=124, file='moment_y'        ,status=file_status,position=file_position,form=file_form)
    open(unit=125, file='moment_z'        ,status=file_status,position=file_position,form=file_form)
    open(unit=223, file='force_x'         ,status=file_status,position=file_position,form=file_form)
    open(unit=224, file='force_y'         ,status=file_status,position=file_position,form=file_form)
    open(unit=225, file='force_z'         ,status=file_status,position=file_position,form=file_form)

    if(.not. lrigid)  then ! integration of vel inside particle wt rigid body assumption
       open(unit=301,file='body_int_u',status=file_status,position=file_position,form=file_form)
       open(unit=302,file='body_int_v',status=file_status,position=file_position,form=file_form)
       open(unit=303,file='body_int_w',status=file_status,position=file_position,form=file_form)
       if(lrotation) then
          open(unit=304,file='body_int_r_u',status=file_status,position=file_position,form=file_form)
          open(unit=305,file='body_int_r_v',status=file_status,position=file_position,form=file_form)
          open(unit=306,file='body_int_r_w',status=file_status,position=file_position,form=file_form)
       endif
    endif
    
    if(ldebug) write(*,*) 'initialization done'
  end subroutine initialization
end module init_m
