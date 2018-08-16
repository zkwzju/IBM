module init_m
  use common_m
  INTEGER ch_fin,m_fin,sph_fin,iters,itfirst,iend,i,j,k
  INTEGER tt_yes,ngbal,n1d
  REAL*8  dt0,dt_p,KE,Ep,diss,ddtEp

contains
  subroutine initialization
    use parser_m

    print *, 'Program running'


    
    !  ch_fin = 0 generates an initial solution
    !  ch_fin = 1 reads the initial data from ch_file.(istart)
    !  outfirst = whether to output data (other than c and u) at first time-step

    call parser_read('ch_file',ch_file)
    call parser_read('sp_file',sp_file)
    call parser_read('rlenx',rlenx)
    call parser_read('rleny',rleny)
    call parser_read('rlenz',rlenz)
    call parser_read('uvbb',uvbb)
    call parser_read('uvbt',uvbt)
    call parser_read('istart',istart)
    call parser_read('iters',iters)
    call parser_read('ch_fin',ch_fin)
    call parser_read('m_fin',m_fin)
    call parser_read('oi_chan',oi_chan)
    call parser_read('out_press',out_press)
    call parser_read('oi_timer',oi_timer)
    call parser_read('oi_cfl',oi_cfl)
    call parser_read('oi_mean',oi_mean)
    call parser_read('oi_spec',oi_spec)
    call parser_read('oi_gbal',oi_gbal)
    call parser_read('oi_1d',oi_1d)
    call parser_read('Reynolds number',re)
    call parser_read('Pressure gradient',ppA)
    call parser_read('Time step dt',dt_p)
    call parser_read('Max CFL',cflmax)
    call parser_read('ppT',ppT)
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

    FLSTLE  = ws.ne.0.0

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
    call initial(dt0,ch_fin,m_fin)

    ! output grid
    open(100,file='grid.dat')
    write(100,*) 'x',nx+1
    write(100,*) xets(1:nx+1)
    write(100,*) 'y',ny+1
    write(100,*) yets(1:ny+1)
    write(100,*) 'z',nz
    write(100,*) zets(1:nz)
    close(100)

    ! open file for velocity difference between fluid and Lagrangian marker
    open(666,file='vel_diff.dat',STATUS='UNKNOWN')
    
    ! produce complete output for initial condition if ...
    if(ch_fin.eq.0)then
       call output

    endif !(ch_fin.eq.0)

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

    open(unit=114, file='liftx')
    open(unit=115, file='drag')
    open(unit=116, file='liftz')
    open(unit=117, file='center_u')
    open(unit=118, file='center_v')
    open(unit=119, file='center_w')
    open(unit=120, file='om_x')
    open(unit=121, file='om_y')
    open(unit=122, file='om_z')
    open(unit=123, file='moment_x')
    open(unit=124, file='moment_y')
    open(unit=125, file='moment_z')

  end subroutine initialization
end module init_m
