      SUBROUTINE advance(code)
      use common_m
      use ellip_common_m
      use p_dyn_m
      implicit none
      
      INTEGER code
c
      REAL*8 cfnl1,cfnl2,cfdif,cfdifpv
      REAL*8 dtcr,idtcr,dtpv,idtpv,prfc1,prfc2,temp,dt0
      INTEGER i,j,jimag,k,irk,nthreads,n
      REAL*8 chbd1t(nz),chbd1b(nz),tmp(nz)
      REAL*8 auxeros(nx0,ny0)
      REAL*8 ddot
      real*8 t_step
      EXTERNAL ddot
      
c
#ifdef OPENMP
      INTEGER  OMP_GET_NUM_THREADS
      EXTERNAL OMP_GET_NUM_THREADS
#endif
c
      if (FLTIMER.and.code.eq.0) call begin_timer(40)
c
      do irk = 1,3
        irkk=irk
        if (FLTIMER) call begin_timer(42)

        if(lcheck .or. (.not. lrigid)) then
           u_temp = u
           v_temp = v
           w_temp = w
           call fft_f2r_2d_new(u_temp)
           call fft_f2r_2d_new(v_temp)
           call fft_f2r_2d_new(w_temp)
        end if
        
        cfnl1   = cfnl1v(irk)
        cfnl2   = cfnl2v(irk)
        cfdif   = cfdifv(irk)
        dtcr    = cfdif*2.d0
        idtcr   = 1.d0/dtcr
c
        cfdifpv = cfdifv(mod(irk+1,3)+1)
c     for restarting with different time step
        if (istart.gt.0.and.itime.eq.istart+1.and.irk.eq.1)
     &   cfdifpv= dt_old/8.d0
        dtpv    = cfdifpv*2.d0
        idtpv   = 1.d0/dtpv !=0.5d0/cfdifpv
c
        prfc1  = 1.d0+(dtcr/dtpv)
        prfc2  = dtcr/dtpv
c
c     compute variable pressure gradient
        if(FLPGVAR)gradpbar=ppA
     &     *dcos(cfpg*(rtime+(cfcumv(irk)-cfdelv(irk))*dt))
c     &     *dsin(cfpg*(rtime+(cfcumv(irk)-cfdelv(irk))*dt))
c
c     compute erosion at bottom now that we have velocity
cNOTE:in this way, erosion bc is put explicitely as afunction of the 
c     previous time shear stress 
c     we may need to compute first the velocities and the compute erosion
c     to do this, there is a problem with the volumetric coupling
c     since div(u)=f(tt) 
c     another idea is to extrapolate the erosion bc
c        if (FLSTLE.and.ttgb(1).eq.-1.d0) then !compute erosion at bottom
c           call comp_erosion(auxeros) !erosion returns in real space
c           call fft_r2f_2d_new_slice(auxeros) !transform to Forier space
c        endif
c
c     compute h of eqns for u, v and w 
        call comp_h(cfnl1,irk)
c
c     compute h of eqn for tt
        if(FLTHRM) call tt_comp_h(cfnl1)
c     compute the right hand sides for u, v and w, results stored in u, v and w
        call rhs(cfnl2,cfdif)
c     compute the right hand side for tt, result stored in tt
        if (FLTHRM) call tt_rhs(cfnl2,cfdif)
c
c     update bc's in Fourier space
c     velocity bc's
C$OMP PARALLEL DEFAULT(SHARED) SHARED(nthreads) 
#ifdef OPENMP
        nthreads = OMP_GET_NUM_THREADS()
#else
        nthreads = 1
#endif
C$OMP DO PRIVATE(i,j,jimag) 
        do j = 1,nyh
          jimag=nyh+j
          do i = 1,nx
c u top bc
           ubctax(i,j    ) = (prfc1*pbxnt(i,j)-prfc2*pbxot(i,j))*dtcr
     &                        + ubct(i,j    )
           ubctax(i,jimag) = (prfc1*pcxnt(i,j)-prfc2*pcxot(i,j))*dtcr
     &                        + ubct(i,jimag)
c u bottom bc
           ubcbax(i,j    ) = (prfc1*pbxnb(i,j)-prfc2*pbxob(i,j))*dtcr
     &                        + ubcb(i,j    )
           ubcbax(i,jimag) = (prfc1*pcxnb(i,j)-prfc2*pcxob(i,j))*dtcr
     &                        + ubcb(i,jimag)
c v top bc
           vbctax(i,j    ) = (prfc1*pbynt(i,j)-prfc2*pbyot(i,j))*dtcr
     &                        + vbct(i,j    )
           vbctax(i,jimag) = (prfc1*pcynt(i,j)-prfc2*pcyot(i,j))*dtcr
     &                        + vbct(i,jimag)
c v bottom bc
           vbcbax(i,j    ) = (prfc1*pbynb(i,j)-prfc2*pbyob(i,j))*dtcr
     &                        + vbcb(i,j    )
           vbcbax(i,jimag) = (prfc1*pcynb(i,j)-prfc2*pcyob(i,j))*dtcr
     &                        + vbcb(i,jimag)
          enddo
        enddo
C$OMP END DO
C$OMP END PARALLEL
c
c     tt bc's
cNOTE:tt bc's - Implemented only erosion at bottom
c     ttab(1)=0.d0, ttbb(1)=1.d0 and ttgb(1)=1.d0
c        if (FLSTLE.and.ttgb(1).eq.-1.d0) then !compute erosion at bottom
cc           call comp_erosion(auxeros) !erosion returns in real space
cc           call fft_r2f_2d_new_slice(auxeros) !transform to Forier space
cC$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,jimag) 
c           do j = 1,nyh
c              jimag=nyh+j
c              do i = 1,nx
cc top bc
cc                ttbct(i,j    )=
cc                ttbct(i,jimag)=
cc bottom bc
cc                ttbcb(i,j    )=-re*pr*abs(w_settling)*auxeros(i,j    )
cc                ttbcb(i,jimag)=-re*pr*abs(w_settling)*auxeros(i,jimag)
c              enddo
c           enddo
cC$OMP END PARALLEL DO
c        endif !FLTHRM.and.ttgb(1).eq.1.d0
c
c     solve the helmholtz equations for u, v, w and tt
        if (mod(itime-istart,oi_cfl).eq.0.and.irk.eq.3) !prnthreads = 1
     &       print *, 'helmholz: time, nthreads = ',itime, nthreads
c     solve the helmholtz equation for u
        call gen_helmholz(re/cfdif,u,ubctax,ubcbax,
     &       uve_mimi,uvci_mimi,uvcn_mimi,
     &       uvprbc1,uvmbct1,uvmbcb1,uvprbcn,uvmbctn,uvmbcbn)
c     solve the helmholtz equation for v
        call gen_helmholz(re/cfdif,v,vbctax,vbcbax,
     &       uve_mimi,uvci_mimi,uvcn_mimi,
     &       uvprbc1,uvmbct1,uvmbcb1,uvprbcn,uvmbctn,uvmbcbn)
c     solve the helmholtz equation for w
        call gen_helmholz(re/cfdif,w,wbct,wbcb,
     &       we_mimi,wci_mimi,wcn_mimi,
     &       wprbc1,wmbct1,wmbcb1,wprbcn,wmbctn,wmbcbn)
c     solve the helmholtz equation for tt
        if (FLTHRM) call gen_helmholz(re*pr/cfdif,tt,ttbct,ttbcb,
     &       tte_mimi,ttci_mimi,ttcn_mimi,
     &       ttprbc1,ttmbct1,ttmbcb1,ttprbcn,ttmbctn,ttmbcbn)

        !>>>
!     particle dynamics before pressure correction for fluid
        if(flagibm  .and. frac_step>0.D0) then
           t_step = frac_step*dtcr
           
           if(ibm_moving .eq.1) call particle_motion(dt,t_step,irk)
           if(lcheck)   call cal_u_diff(irk)
           if(.not. lrigid) then
              do n=1,num_p
                 call part_inner_integral(u_temp,v_temp,w_temp,
     &                int_u1(n),int_v1(n),int_w1(n),
     &                int_r_u1(n),int_r_v1(n),int_r_w1(n),n)
              enddo
              call body_integration_output(irk)
           endif
        endif
        !<<<
        
c     compute the divergence for the pressure step
c     the right hand side of eqn for p is stored in rw
        call ddz(w,rw)!ddz computes in real and Fourier spaces
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,jimag,k,temp),
C$OMP$                             SHARED(idtcr)
        do k=1,nz
          do j=1,nyh
            jimag=nyh+j
c add dw/dz
            call dscal(nx,idtcr,rw(1,j    ,k),1)
            call dscal(nx,idtcr,rw(1,jimag,k),1)
c add dv/dy
            temp = idtcr*wavey(j)
            call daxpy(nx,-temp,v(1,jimag,k),1,rw(1,j,    k),1)
            call daxpy(nx, temp,v(1,j    ,k),1,rw(1,jimag,k),1)
          enddo
c add du/dx
          do i=1,nx
            temp = idtcr*wavex(i)
            call daxpy(nyh,-temp,u(i,nyhp,k),nx0,rw(i,1   ,k),nx0)
            call daxpy(nyh, temp,u(i,1   ,k),nx0,rw(i,nyhp,k),nx0)
          enddo
        enddo
C$OMP END PARALLEL DO
c
c     solve the helmholtz equation for pressure p
c     rhs for pressure equation is in rw
c     pressure returned in Fourie space in variable p
        call pstep
c
c     compute u(i) from u(*i) and dp/dx
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,k,temp)
        do k=1,nz
          do i=1,nx
            temp = dtcr*wavex(i)
            call daxpy(nyh, temp,p(i,nyhp,k),nx0, u(i,1,k),nx0)
            call daxpy(nyh,-temp,p(i,1,k),nx0,u(i,nyhp,k),nx0)
          enddo
        enddo
C$OMP END PARALLEL DO
c
c     compute v(i) from v(*i) and dp/dy
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(j,jimag,k,temp)
        do k=1,nz
          do j=1,nyh
            jimag=nyh+j
c
            temp = dtcr*wavey(j)
            call daxpy(nx, temp,p(1,jimag,k),1, v(1,j,k),1)
            call daxpy(nx,-temp,p(1,j,k),1,v(1,jimag,k),1)
          enddo
        enddo
C$OMP END PARALLEL DO
c
c     compute w(i) from w(*i) and dp/dz
        call ddz(p,ru)!ddz computes in real and Fourier spaces
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(j,jimag,k)
        do k = 1,nz
          do j=1,nyh
            jimag=nyh+j
            call daxpy(nx,-dtcr,ru(1,j,k),1, w(1,j,k),1)
            call daxpy(nx,-dtcr,ru(1,jimag,k),1,w(1,jimag,k),1)
          enddo
cjpf  zero out the mean of w (which can grow very slowly over time)
cc          w( 1,1,k) = 0.d0
cc          w(1,nyhp,k) = 0.d0
        enddo
C$OMP END PARALLEL DO
c
        if (FLTIMER) call end_timer(42)
c
        if (FLTIMER) call begin_timer(43)
c     update pressure gradient for bc's for tangencial velocity
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,jimag)
        do j = 1,nyh
           jimag=nyh+j
           do i = 1,nx
c top 
              pbxot(i,j) = pbxnt(i,j)
              pcxot(i,j) = pcxnt(i,j)
              pbyot(i,j) = pbynt(i,j)
              pcyot(i,j) = pcynt(i,j)
c bottom
              pbxob(i,j) = pbxnb(i,j)
              pcxob(i,j) = pcxnb(i,j)
              pbyob(i,j) = pbynb(i,j)
              pcyob(i,j) = pcynb(i,j)
           enddo
        enddo
C$OMP END PARALLEL DO
c
c update values of aux pressure at boundaries for aux vel bc's
c if needed ddz is computed, but ddx and ddy are computed below
c top boundary 
        if(uvat.eq.1.d0)then
c no-slip
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(j,jimag)
           do j = 1,nyh
              jimag=nyh+j
              call dcopy(nx,p(1,jimag,1),1,pbxnt(1,j),1)
              call dcopy(nx,p(1,jimag,1),1,pbynt(1,j),1)
              call dcopy(nx,p(1,j    ,1),1,pcxnt(1,j),1)
              call dcopy(nx,p(1,j    ,1),1,pcynt(1,j),1)
           enddo
C$OMP END PARALLEL DO
        elseif(uvbt.eq.1.d0)then
c slip
           call dcopy(nz,chbd1( 1,1),nz0,chbd1t,1)
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,jimag,tmp)
           do j = 1,nyh
              jimag=nyh+j
              do i=1,nx
c     ddz and place into imag part
                 call dcopy(nz,p(i,j,1),nx0y0,tmp,1)
                 pcynt(i,j) = ddot(nz,chbd1t,1,tmp,1)
                 pcxnt(i,j) = pcynt(i,j)
c     ddz and place into real part
                 call dcopy(nz,p(i,jimag,1),nx0y0,tmp,1)
                 pbynt(i,j) = ddot(nz,chbd1t,1,tmp,1)
                 pbxnt(i,j) = pbynt(i,j)
              enddo
           enddo
C$OMP END PARALLEL DO
        else
           write(26,*)'Wrong BC for velocity for top boundary'
           write(26,*)'uvbt:',uvbt,'uvat:',uvat
           write( *,*)'Wrong BC for velocity for top boundary'
           write( *,*)'uvbt:',uvbt,'uvat:',uvat
           stop
        endif
c
c bottom boundary
        if(uvab.eq.1.d0)then
c no-slip
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(j,jimag)
           do j = 1,nyh
              jimag=nyh+j
              call dcopy(nx,p(1,jimag,nz),1,pbxnb(1,j),1)
              call dcopy(nx,p(1,jimag,nz),1,pbynb(1,j),1)
              call dcopy(nx,p(1,j    ,nz),1,pcxnb(1,j),1)
              call dcopy(nx,p(1,j    ,nz),1,pcynb(1,j),1)
           enddo
C$OMP END PARALLEL DO
        elseif(uvbb.eq.1.d0)then
c slip
           call dcopy(nz,chbd1(nz,1),nz0,chbd1b,1)
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,jimag,tmp)
           do j = 1,nyh
              jimag=nyh+j
              do i=1,nx
c     ddz and place into imag part
                 call dcopy(nz,p(i,j,1),nx0y0,tmp,1)
                 pcynb(i,j) = ddot(nz,chbd1b,1,tmp,1)
                 pcxnb(i,j) = pcynb(i,j)
c     ddz and place into real part
                 call dcopy(nz,p(i,jimag,1),nx0y0,tmp,1)
                 pbynb(i,j) = ddot(nz,chbd1b,1,tmp,1)   
                 pbxnb(i,j) = pbynb(i,j)
              enddo
           enddo
C$OMP END PARALLEL DO
        else
           write(26,*)'Wrong BC for velocity for bottom boundary'
           write(26,*)'uvbb:',uvbb,'uvab:',uvab
           write( *,*)'Wrong BC for velocity for bottom boundary'
           write( *,*)'uvbb:',uvbb,'uvab:',uvab
           stop
        endif
c
c     Compute ddx and ddy of aux pressure for aux vel bc's
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,jimag)
        do j = 1,nyh
           jimag=nyh+j
           do i = 1,nx
c top
              pbxnt(i,j) = -wavex(i)*pbxnt(i,j)
              pcxnt(i,j) =  wavex(i)*pcxnt(i,j)
              pbynt(i,j) = -wavey(j)*pbynt(i,j)
              pcynt(i,j) =  wavey(j)*pcynt(i,j)
c bottom
              pbxnb(i,j) = -wavex(i)*pbxnb(i,j)
              pcxnb(i,j) =  wavex(i)*pcxnb(i,j)
              pbynb(i,j) = -wavey(j)*pbynb(i,j)
              pcynb(i,j) =  wavey(j)*pcynb(i,j)
           enddo
        enddo
C$OMP END PARALLEL DO
        if (FLTIMER) call end_timer(43)
c     
c     compute pressure (needed for transport of TKE)
c        if(out_press.ne.0.and.
c     &       mod((itime-istart),oi_chan).eq.0.and.irk.eq.3)then
c           call viscxyz_f(p,rw,re,FLSV_NO)
cC$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k)
c           do k=1,nz
c              do j=1,ny
c                 do i=1,nx
c                    p(i,j,k)=p(i,j,k)-cfdif*rw(i,j,k)
cc                    pbefore(i,j,k)=p(i,j,k)
c                 enddo
c              enddo
c           enddo
cC$OMP END PARALLEL DO


        if(irk.eq.3)then
           call viscxyz_f(p,rw,re,FLSV_NO)
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k)
           do k=1,nz
              do j=1,ny
                 do i=1,nx
                    p(i,j,k)=p(i,j,k)-cfdif*rw(i,j,k)
                 enddo
              enddo
           enddo
C$OMP END PARALLEL DO          
        endif

!>>>  
!     particle dynamics after pressure correction for fluid
        if(flagibm  .and. frac_step<1.0D0) then
           t_step = (1.0D0-frac_step)*dtcr
           
           if(ibm_moving .eq.1)  call particle_motion(dt,t_step,irk)

           if(lcheck)   call cal_u_diff(irk)
           if(.not. lrigid) then
              do n=1,num_p
                 call part_inner_integral(u_temp,v_temp,w_temp,
     &                int_u1(n),int_v1(n),int_w1(n),
     &                int_r_u1(n),int_r_v1(n),int_r_w1(n),n)
                 call body_integration_output(irk)
              enddo
            endif
        endif
!<<<
      enddo                     !irk
c     
      if (FLTIMER) call end_timer(40)
c     
      return
      END

! -------------------------------------------------- !
! -------------------------------------------------- !
      subroutine cal_u_diff(irk)
      use common_m
      use ellipsoid_m
      use ellip_common_m

      implicit none

      integer irk
      integer l,n,i1,i2,j1,j2,k1,k2
      integer i,j,k,i_t,j_t,i_ddf,j_ddf
      integer nblock,nremain,nn
      

      if(irk .eq. 1 .or. irk .eq. 2) return  ! write if only irk=3

      diff_v = 0.0D0
      
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(l,n,i1,i2,j1,j2,k1,k2
C$OMP& ,i,j,k,i_t,j_t,i_ddf,j_ddf)
      do n=1,num_p
         do l=1,n_l(n)
            i1=p_iw(l,n)
            i2=p_ie(l,n)
            j1=p_js(l,n)
            j2=p_jn(l,n)
            k1=p_kb(l,n)
            k2=p_kt(l,n)
            do i=i1,i2
               i_ddf=i-i1+1
               if(i .lt. 1) then
                  i_t=nx+i
               elseif (i .gt. nx) then
                  i_t=i-nx
               else
                  i_t=i
               endif
               do j=j1,j2
                  j_ddf=j-j1+1
                  if(j .lt.1) then
                     j_t=ny+j
                  elseif(j .gt. ny) then
                     j_t=j-ny
                  else
                     j_t=j
                  endif
                  do k=k1,k2
                     diff_v(l,n)=diff_v(l,n)+v_temp(i_t,j_t,k)
     &                    *ddf(l,n,i_ddf,j_ddf,k-k1+1)
                  enddo
               enddo
            enddo
            
            diff_v(l,n)=diff_v(l,n)-v_p(l,n) ! vel difference

         enddo
      enddo
C$OMP END PARALLEL DO
         
      if(lunformatted) then
         write(666) rtime,diff_v
      else
         nblock=ceiling(n_ll*1.d0/MaxPerLine)
         nremain=mod(n_ll,MaxPerLine)
         do n=1,num_p
            write(666,'(ES16.8)',advance='no') rtime
            do nn=1,nblock
               if(nn .eq. nblock .and. nremain .ne. 0) then
                  write(666,'(<nremain>ES16.8)') (diff_v(l,n),
     &                 l=(nn-1)*MaxPerLine+1,n_l(n))
               else
                  write(666,'(<MaxPerLine>ES16.8)') (diff_v(l,n),
     &                 l=(nn-1)*MaxPerLine+1,nn*MaxPerLine)
               endif
            enddo
         enddo
      endif
      
      end subroutine

! -------------------------------------------------- !
      subroutine body_integration_output(irk)
      use common_m
      use ellipsoid_m
      use ellip_common_m
      implicit none

      integer irk
      integer nblock,nremain,nn,l

      if(irk .eq. 1 .or. irk .eq. 2) return  ! write if only irk=3
      nblock=ceiling(num_p/100.0d0)
      nremain=mod(num_p,100)
      write(301,'(I10)',advance='no') itime
      write(302,'(I10)',advance='no') itime
      write(303,'(I10)',advance='no') itime

      do nn=1,nblock
         if(nn .eq. nblock) then
            write(301,'(<nremain>ES16.8)')  (int_u1(l),
     &           l=(nn-1)*100+1,num_p)
            write(302,'(<nremain>ES16.8)')  (int_v1(l),
     &           l=(nn-1)*100+1,num_p)
            write(303,'(<nremain>ES16.8)')  (int_w1(l),
     &           l=(nn-1)*100+1,num_p)
         else
            write(301,'(100ES16.8)') (int_u1(l),
     &           l=(nn-1)*100+1,nn*100)
            write(302,'(100ES16.8)') (int_v1(l),
     &           l=(nn-1)*100+1,nn*100)
            write(303,'(100ES16.8)') (int_w1(l),
     &           l=(nn-1)*100+1,nn*100)
         endif
      enddo
      if(lrotation) then
         write(304,'(I10)') itime
         write(305,'(I10)') itime
            write(306,'(I10)') itime
            do nn=1,nblock
               if(nn .eq. nblock) then
                  write(304,'(<nremain>ES16.8)')  (int_r_u1(l),
     &                 l=(nn-1)*100+1,num_p)
                  write(305,'(<nremain>ES16.8)')  (int_r_v1(l),
     &                 l=(nn-1)*100+1,num_p)
                  write(306,'(<nremain>ES16.8)')  (int_r_w1(l),
     &                 l=(nn-1)*100+1,num_p)
               else
                  write(304,'(100ES16.8)') (int_r_u1(l),
     &                 l=(nn-1)*100+1,nn*100)
                  write(305,'(100ES16.8)') (int_r_v1(l),
     &                 l=(nn-1)*100+1,nn*100)
                  write(306,'(100ES16.8)') (int_r_w1(l),
     &                 l=(nn-1)*100+1,nn*100)
               endif
            enddo
         endif

      end subroutine body_integration_output
! -------------------------------------------------- !
      subroutine weighted_fluid_vel(weight_v)
      use common_m
      use ellip_common_m,only:x_0,x_1,y_0,y_1,z_0,z_1
      implicit none
      real*8 weight_v
      integer n,l,i,j,k,id,jd,i1,i2,j1,j2,k1,k2

      weight_v = 0.d0
      
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(n,i1,i2,j1,j2,k1,k2,l,i)
!$OMP$  private(id,j,jd,k,weight_v)
      do n=1,num_p
         i1 = x_0(n)
         i2 = x_1(n)
         j1 = y_0(n)
         j2 = y_1(n)
         k1 = z_0(n)
         k2 = z_1(n)

         do l=1,n_l(n)
            do i=i1,i2
               if(i<1) then
                  id = i+nx
               elseif(i>nx) then
                  id = i-nx
               else
                  id = i
               end if
               
               do j=j1,j2
                  if(j<1) then
                     jd = j+ny
                  elseif(i>nx) then
                     jd = j-ny
                  else
                     jd = j
                  end if
                  do k=k1,k2
                     weight_v=weight_v+
     &                    ddf(l,n,i-i1+1,j-j1+1,k-k1+1)*
     &                    cell_v(k)*v_temp(id,jd,k)
                  end do 
               end do 
            end do
         end do
      end do 
!$OMP END PARALLEL DO    
      end subroutine weighted_fluid_vel
      

! --------------------------------------------------!
      subroutine particle_motion(delta_t,t_step,irk)
      use common_m
      use ellipsoid_m
      use ellip_common_m
      use p_dyn_m
      implicit none
      
      real*8 delta_t,t_step
      real*8 F_v(3),torq_v(3),it_step
      real*8 xp3(3),up3(3),ax1(3),ax2(3),ax3(3)
      real*8 weight_v
      integer irk,n,l,i,j,k,i_ind,j_ind,k_ind
      real*8 v_diff(MaxPerLine)     ! hardcoded matrix dimension
      
      it_step = 1.0D0/t_step
      
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(n,om_ellip,ax1,ax2,ax3,up3),
!$OMP& private(I_ellip,F_v,torq_v,om_ellip_b)
      do n=1,num_p

         om_ellip(1) = om_x(n)
         om_ellip(2) = om_y(n)
         om_ellip(3) = om_z(n)
         
         ax1 = axis_1(:,n)
         ax2 = axis_2(:,n)
         ax3 = axis_3(:,n)

         xp3(1) = x_c(n)
         xp3(2) = y_c(n)
         xp3(3) = z_c(n)
            
         up3(1) = u_c(n)
         up3(2) = v_c(n)
         up3(3) = w_c(n)

         I_ellip = I_particle(:,n)

         if(lrigid) then
           
            if(ltranslation) then
               F_v(1) = - for_px(n)*(delta_t/t_step) ! stupid definition
               F_v(2) = - for_py(n)*(delta_t/t_step) ! in direct forcing 
               F_v(3) = - for_pz(n)*(delta_t/t_step)

               if(lclip_x) then
                  F_v(1) = 0.0D0
               endif
               if(lclip_y) then
                  F_v(2) = 0.0D0
                  grav   = 0.0D0
               endif
               if(lclip_z) then
                  F_v(3) = 0.0D0
               endif

               call part_translation2(xp3,up3,F_v,grav,vol_ellip,t_step)
               
            endif
            if(lrotation) then
               torq_v(1) = - torq_x(n)*(delta_t/t_step) ! stupid definition
               torq_v(2) = - torq_y(n)*(delta_t/t_step) ! in direct forcing
               torq_v(3) = - torq_z(n)*(delta_t/t_step) ! 
               if(lclip_ox) then
                  torq_v(1) = 0.0D0
               endif
               if(lclip_oy) then
                  torq_v(2) = 0.0D0
               endif
               if(lclip_oz) then
                  torq_v(3) = 0.0D0
               endif
               if(lsphere) then
                  call part_rotation2(om_ellip, torq_v,I_ellip,t_step)
               else
                  torq_v = rho_f*rho_p/(rho_p-rho_f)*torq_v
                  call part_rotation(ax1,ax2,ax3,om_ellip,om_ellip_b,
     &                 torq_v,I_ellip,t_step)
               endif
            endif
         else
            if(ltranslation) then
               F_v(1) = (int_u1(n)-int_u0(n))*it_step
               F_v(2) = (int_v1(n)-int_v0(n))*it_step
               F_v(3) = (int_w1(n)-int_w0(n))*it_step

               F_v(1) = F_v(1) - for_px(n)*(delta_t/t_step)
               F_v(2) = F_v(2) - for_py(n)*(delta_t/t_step)
               F_v(3) = F_v(3) - for_pz(n)*(delta_t/t_step)

               if(lclip_x) then
                  F_v(1) = 0.0D0
               endif
               if(lclip_y) then
                  F_v(2) = 0.0D0
                  grav   = 0.0D0
               endif
               if(lclip_z) then
                  F_v(3) = 0.0D0
               endif
               call part_translation(xp3,up3,F_v,grav,vol_ellip,t_step)
            endif

            if(lrotation) then
               torq_v(1) = (int_r_u1(n)-int_r_u0(n))*it_step
               torq_v(2) = (int_r_v1(n)-int_r_v0(n))*it_step
               torq_v(3) = (int_r_w1(n)-int_r_w0(n))*it_step

               torq_v(1) = torq_v(1) - torq_x(n)*(delta_t/t_step)
               torq_v(2) = torq_v(2) - torq_y(n)*(delta_t/t_step)
               torq_v(3) = torq_v(3) - torq_z(n)*(delta_t/t_step)
               if(lclip_ox) then
                  torq_v(1) = 0.0D0
               endif
               if(lclip_oy) then
                  torq_v(2) = 0.0D0
               endif
               if(lclip_oz) then
                  torq_v(3) = 0.0D0
               endif
               if(lsphere) then
                  call part_rotation2(om_ellip, torq_v,I_ellip,t_step)
               else
                  call part_rotation(ax1,ax2,ax3,om_ellip,om_ellip_b,
     &                 torq_v,I_ellip,t_step)
               endif
            endif
            int_u0(n) = int_u1(n)
            int_v0(n) = int_v1(n)
            int_w0(n) = int_w1(n)
            int_r_u0(n) = int_r_u1(n)
            int_r_v0(n) = int_r_v1(n)
            int_r_w0(n) = int_r_w1(n)
         endif
               
  
         x_c(n) = xp3(1)
         y_c(n) = xp3(2)
         z_c(n) = xp3(3)
         
         u_c(n) = up3(1)
         v_c(n) = up3(2)
         w_c(n) = up3(3)
         

         axis_1(:,n) = ax1
         axis_2(:,n) = ax2
         axis_3(:,n) = ax3
         
         om_x(n) = om_ellip(1)
         om_y(n) = om_ellip(2)
         om_z(n) = om_ellip(3)
         
      enddo
!$OMP END PARALLEL DO  
      return
      end subroutine particle_motion
