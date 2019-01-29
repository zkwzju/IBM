module global_m
  implicit none
  INTEGER nx,ny,nz ! Number of grid points in each direction
  INTEGER num_p,n_ll ! 

  INTEGER np,dnp !recall that if np=1 => dnp=1, number of phases to average
  INTEGER m ! Number of scalar fields

  !  Derived parameters
  INTEGER nxl,nyl
  INTEGER nxh,nyh,nxhp,nyhp,nxhm,cnxh,nzm,nzmm
  INTEGER nxlh,nylh,nxlhp,nylhp
  INTEGER nxlch,nxlchm,nxlchp,cnxlch,cnxlchm
  INTEGER nylch,nylchp

  !  Non-final array subscript sizes (should avoid powers of two)
  INTEGER nx0,ny0,nz0,nzmm0,nxhp0,nyhp0,tnxhp0,tnyhp0
  INTEGER nxl0,nyl0,nxlhp0,nylhp0,tnxlhp0,tnylhp0
  INTEGER nxny,nx0y,nx0ylch,nx0y0,nx0y0z,nz0z
  INTEGER nyhp0x0,nylhp0xl0,nxl0yl

contains

  subroutine global_init
    use parser_m
    integer :: iunit=345,istatus

    character(len=120) :: input,command_buffer
    external getarg
    ! Read inputs from command line
    if (iargc()==1) then
       ! here input file defaults
       call getarg(1,command_buffer)
       read(command_buffer,*) input 
    else
       input='LP_input'
    end if

    call parser_init
    call parser_parsefile(input,iunit, istatus)
    if(istatus .eq. 1) then    
       write(*,*) 'Fatal: Settings file LP_input not found'
       stop
    end if

    call parser_read('Grid nx',nx,100)
    call parser_read('Grid ny',ny,100)
    call parser_read('Grid nz',nz,100)
    call parser_read('Number of particle',num_p,1)
    call parser_read('Max Lagrangian markers',n_ll,1000)

    np=1
    dnp=1

    m=1
    nxl=3*nx/2
    nyl=3*ny/2
    nxh = nx/2
    nyh = ny/2
    nxhp = nxh+1
    nyhp = nyh+1
    nxhm = nxh-1
    cnxh = nxh+2
    nxlh = nxl/2
    nylh = nyl/2
    nxlhp = nxlh+1
    nylhp = nylh+1
    nzm = nz-1
    nzmm=nz-2
    nxlch = nxh
    nxlchm = nxlch-1
    nxlchp = nxlch+1
    cnxlch = nxl+2-nxlch
    cnxlchm = cnxlch-1
    nylch = nyh
    nylchp = nylch+1
    nx0=nx+1
    ny0=ny+1
    nz0=nz+1
    nzmm0=nzmm+1
    nxhp0=nxhp
    nyhp0=nyhp
    tnxhp0=2*nxhp0
    tnyhp0=2*nyhp0
    nxl0=nxl+1
    nyl0=nyl+1
    nxlhp0=nxlhp
    nylhp0=nylhp
    tnxlhp0=2*nxlhp0
    tnylhp0=2*nylhp0
    nxny=nx*ny
    nx0y=nx0*ny
    nx0ylch=nx0*nylch
    nx0y0=nx0*ny0
    nx0y0z=nx0*ny0*nz
    nz0z=nz0*nz
    nyhp0x0=nyhp0*nx0
    nylhp0xl0=nylhp0*nxl0
    nxl0yl=nxl0*nyl

!    write(*,*) 'global_init done'

  end subroutine global_init

end module global_m

