program pv
  use globals
  use pv3D
  use cema, only: calc_cema, allocation_cema, deallocation_cema, EI, &
                  read_species_names_cema, read_indices_cema, species_names_cema
  implicit none
  double precision::dYfdx,dYfdy,dYodx,dYody
!   integer i_tmp
  integer :: ptr_states, ptr_Y, ptr_analyze, ptr_EI
  
  external :: allocation, Z_ini, chem_ini, calc_z, calc_Dh, deallocation, pv_particle_output

  open(200,file='./dat/setting_pv.dat',form='formatted')
   read(200,'()') !--- step infomation ---!
   read(200,*) step0
   read(200,*) step1
   read(200,*) step2
   read(200,'()') !--- domain information ---!
   read(200,*) x_sta
   read(200,*) x_end
   read(200,*) y_sta
   read(200,*) y_end
   read(200,*) z_sta
   read(200,*) z_end
   read(200,*) nx 
   read(200,*) ny
   read(200,*) nz
   read(200,*) ibd
   read(200,*) jbd
   read(200,*) kbd
   read(200,*) iprocs
   read(200,*) jprocs
   read(200,*) kprocs
   read(200,'()') !--- reaction information ---!
   read(200,*) chem_dir
   read(200,*) nf
   read(200,*) nrf
   read(200,*) nrb
   read(200,*) nrp
   read(200,'()') !--- flags ---!
   read(200,*) flag_particle
   read(200,*) flag_Z
   read(200,*) flag_pv    
   read(200,*) flag_Y    
   read(200,*) flag_EI    
   close(200)

   if(z_sta==z_end)then
      flag_plane=1
   else if(x_sta==x_end)then
      flag_plane=2
   else if(y_sta==y_end)then
      flag_plane=3
   else 
      flag_plane=4
   endif

   write(*,*) "step    ",step0, step1, step2
   write(*,*) "x(vis)  ",x_sta,x_end 
   write(*,*) "y(vis)  ",y_sta,y_end 
   write(*,*) "z(vis)  ",z_sta,z_end
   
   write(*,*) "domain  ",nx, ny, nz
   write(*,*) "bd      ",ibd,jbd,kbd  
   write(*,*) "process ",iprocs,jprocs,kprocs
   write(*,*) "nf      ", nf
   write(*,*) "particle",flag_particle
   write(*,*) "plane   ",flag_plane

   ! calculate pointer for sf
   ptr_states = 3
   ptr_Y = ptr_states + 4

   ! whether output Y
   if (flag_Y) then
      ptr_analyze = ptr_Y + nf
   else
      ptr_analyze = ptr_Y
   end if

   ptr_EI = ptr_analyze + 12
   
   ! whether output EI
   if (flag_EI) then
      ns = ptr_EI + nf
   else
      ns = ptr_EI
   end if

  call allocation
  call allocation_cema()
  call read_species_names_cema()
  call read_indices_cema()

  !  allocate(wg(1:nf))
 

  open(20,file='../../xg.dat')                      
  do i=1-ibd,nx+ibd
     read(20,*) ii,xu(i),xg(i),dxu(i),dxg(i)
  end do
  close(20)

  open(20,file='../../yg.dat')
  do j=1-jbd,ny+jbd
    read(20,*) jj,yv(j),yg(j),dyv(j),dyg(j)
  end do
  close(20)

  open(20,file='../../zg.dat')
  do k=1-kbd,nz+kbd
    read(20,*) kk,zw(k),zg(k),dzw(k),dzg(k)
  end do
  close(20)

  if(flag_Z) call Z_ini
  if(flag_Z) call chem_ini
  

  mainloop:do step=step0,step1,step2
     write(filenumber,'(i8.8)')step

     do myrank=0,(iprocs*jprocs*kprocs-1)
        write(cpunumber,'(i5.5)')myrank
        
        m=0
        do k=0,kprocs-1
        do j=0,jprocs-1
        do i=0,iprocs-1
           if(myrank==m)then
              myrank_i=i
              myrank_j=j
              myrank_k=k
           end if
           m=m+1
        end do
        end do
        end do
      
        ista=nx/iprocs*myrank_i+1
        iend=nx/iprocs*(myrank_i+1)
        jsta=ny/jprocs*myrank_j+1
        jend=ny/jprocs*(myrank_j+1)
        ksta=nz/kprocs*myrank_k+1
        kend=nz/kprocs*(myrank_k+1)

        if( .not.(iend<x_sta .or. ista>x_end).and. .not.(jend<y_sta .or. jsta>y_end).and. .not.(kend<z_sta .or. ksta>z_end) )then
           allocate(u_local(ista-ibd:iend+ibd,jsta-jbd:jend+jbd,ksta-kbd:kend+kbd))
           allocate(v_local(ista-ibd:iend+ibd,jsta-jbd:jend+jbd,ksta-kbd:kend+kbd))
           allocate(w_local(ista-ibd:iend+ibd,jsta-jbd:jend+jbd,ksta-kbd:kend+kbd))
           allocate(r_local(ista-ibd:iend+ibd,jsta-jbd:jend+jbd,ksta-kbd:kend+kbd))
           allocate(p_local(ista-ibd:iend+ibd,jsta-jbd:jend+jbd,ksta-kbd:kend+kbd))
           allocate(t_local(ista-ibd:iend+ibd,jsta-jbd:jend+jbd,ksta-kbd:kend+kbd))
           allocate(h_local(ista-ibd:iend+ibd,jsta-jbd:jend+jbd,ksta-kbd:kend+kbd))
           allocate(y_local(ista-ibd:iend+ibd,jsta-jbd:jend+jbd,ksta-kbd:kend+kbd,1:nf))

           ifn=90
           if(iprocs*jprocs*kprocs/=1)then
              open(ifn,file='../'//cpunumber//'/f'//filenumber//'.dat',form='unformatted') !,convert='big_endia)
              write(*,*)'now reading ./'//cpunumber//'/f'//filenumber//'.dat'
           else
              open(ifn,file='../f'//filenumber//'.dat',form='unformatted') !,convert='big_endian') 
              write(*,*)'now reading ./f'//filenumber//'.dat'
           endif

           read(ifn) u_local
           read(ifn) v_local
           read(ifn) w_local
           read(ifn) r_local
           read(ifn) p_local
           read(ifn) t_local
           read(ifn) h_local
           read(ifn) y_local
           close(ifn)
           ! do k=z_sta,z_end !output_all
           ! do j=y_sta,y_end !output_all
           ! do i=x_sta,x_end !output_all

           do k=max(ksta,z_sta),min(kend,z_end)
           do j=max(jsta,y_sta),min(jend,y_end)
           do i=max(ista,x_sta),min(iend,x_end) 
              ! do k=max(ksta,z_sta),min(kend,z_end)
              ! do j=max(jsta,y_sta),min(jend,y_end)
              ! do i=max(ista,x_sta),min(iend,x_end) 
              
              sf(i,j,k,1)=u_local(i,j,k)
              sf(i,j,k,2)=v_local(i,j,k)
              sf(i,j,k,3)=w_local(i,j,k)

              ! ptr_states
              sf(i,j,k,ptr_states+1)=r_local(i,j,k)
              sf(i,j,k,ptr_states+2)=p_local(i,j,k)
              sf(i,j,k,ptr_states+3)=t_local(i,j,k)
              sf(i,j,k,ptr_states+4)=h_local(i,j,k)

              ! ptr_Y
              if (flag_Y) then
               ! asign Y to sf
               do kk = 1, nf
                sf(i,j,k,ptr_Y+kk)=y_local(i,j,k,kk)
               end do
              end if

              ! ptr_analyze
              sf(i,j,k,ptr_analyze+7)=sum(y_local(i,j,k,:))     !sumY

              dyodx=(y_local(i+1,j,k, 3)-y_local(i-1,j,k, 3))/(dxg(i)*dxg(i+1))
              dyody=(y_local(i,j+1,k, 3)-y_local(i,j-1,k, 3))/(dyg(j)*dyg(j+1))
              dyfdx=(y_local(i+1,j,k,11)-y_local(i-1,j,k,11))/(dxg(i)*dxg(i+1))
              dyfdy=(y_local(i,j+1,k,11)-y_local(i,j-1,k,11))/(dyg(j)*dyg(j+1))
              sf(i,j,k,ptr_analyze+8)=(dyodx*dyfdx+dyody*dyfdy+1d-12) &
                                     /(sqrt(dyodx*dyodx+dyody*dyody)*sqrt(dyfdx*dyfdx+dyfdy*dyfdy)+1d-12)
              
              call calc_z(y_local(i,j,k,:),sf(i,j,k,ptr_analyze+1),sf(i,j,k,ptr_analyze+3),&
                                           sf(i,j,k,ptr_analyze+4),sf(i,j,k,ptr_analyze+5))
              call calc_Dh(y_local(i,j,k,:),r_local(i,j,k),t_local(i,j,k),sf(i,j,k,ptr_analyze+6))
              call calc_cema(y_local(i,j,k,:),p_local(i,j,k),t_local(i,j,k),sf(i,j,k,ptr_analyze+9), &
                             sf(i,j,k,ptr_analyze+10), sf(i,j,k,ptr_analyze+11), sf(i,j,k,ptr_analyze+12))

              ! ptr_EI
              if (flag_EI) then
               ! asign EI to sf
               do kk = 1, nf
                  sf(i,j,k,ptr_EI + kk) = EI(kk)
               end do
              end if

           end do
           end do
           end do
           
           call deallocation
         !   call deallocation_cema
        end if
        
     end do

     sf(:,:,:,ptr_analyze+2)=0d0

     if(nz==1)then
        k=z_sta
        do j=y_sta+1,y_end-1
           do i=x_sta+1,x_end-1
              dumdx=(dxg(i)*dxg(i)*sf(i+1,j,k,ptr_analyze+1) &
                   + (dxg(i+1)*dxg(i+1)-dxg(i)*dxg(i))*sf(i,j,k,ptr_analyze+1) &
                   - dxg(i+1)*dxg(i+1)*sf(i-1,j,k,ptr_analyze+1)) &
                   / (dxg(i)*dxg(i+1)*(dxg(i)+dxg(i+1))) 
              dumdy=(dyg(j)*dyg(j)*sf(i,j+1,k,ptr_analyze+1) &
                   + (dyg(j+1)*dyg(j+1)-dyg(j)*dyg(j))*sf(i,j,k,ptr_analyze+1) &
                   - dyg(j+1)*dyg(j+1)*sf(i,j-1,k,ptr_analyze+1)) &
                   / (dyg(j)*dyg(j+1)*(dyg(j)+dyg(j+1))) 
              sf(i,j,k,ptr_analyze+2)=2d0 * sf(i,j,k,ptr_analyze+6) * (dumdx*dumdx+dumdy*dumdy)
           enddo
        enddo
     else
        do k=z_sta+1,z_end-1
           do j=y_sta+1,y_end-1
              do i=x_sta+1,x_end-1
                 dumdx=(dxg(i)*dxg(i)*sf(i+1,j,k,ptr_analyze+1) &
                      + (dxg(i+1)*dxg(i+1)-dxg(i)*dxg(i))*sf(i,j,k,ptr_analyze+1) &
                      - dxg(i+1)*dxg(i+1)*sf(i-1,j,k,ptr_analyze+1)) &
                      / (dxg(i)*dxg(i+1)*(dxg(i)+dxg(i+1))) 
                 dumdy=(dyg(j)*dyg(j)*sf(i,j+1,k,ptr_analyze+1) &
                      + (dyg(j+1)*dyg(j+1)-dyg(j)*dyg(j))*sf(i,j,k,ptr_analyze+1) &
                      - dyg(j+1)*dyg(j+1)*sf(i,j-1,k,ptr_analyze+1)) &
                      / (dyg(j)*dyg(j+1)*(dyg(j)+dyg(j+1))) 
                 dumdz=(dzg(k)*dzg(k)*sf(i,j,k+1,ptr_analyze+1) &
                      + (dzg(k+1)*dzg(k+1)-dzg(k)*dzg(k))*sf(i,j,k,ptr_analyze+1) &
                      - dzg(k+1)*dzg(k+1)*sf(i,j,k-1,ptr_analyze+1)) &
                      / (dzg(k)*dzg(k+1)*(dzg(k)+dzg(k+1))) 
                 sf(i,j,k,ptr_analyze+2)=2d0 * sf(i,j,k,ptr_analyze+6) * (dumdx*dumdx+dumdy*dumdy+dumdz*dumdz)
              enddo
           enddo
        enddo
     endif


     open(333,file='./dat/Zst_ypos_'//filenumber//'.dat',form='formatted')
     open(444,file='./dat/Zst_yneg_'//filenumber//'.dat',form='formatted')
     if(nz==1)then
     write(333,*) " x[m] y[m] chi[1/s] u[m/s] v[m/s] T[K] tau_u[s] tau_mix[s]"
     write(444,*) " x[m] y[m] chi[1/s] u[m/s] v[m/s] T[K] tau_u[s] tau_mix[s]"
     k=z_sta
     tau_up=0d0
     tau_un=0d0
        do i=x_sta,x_end
           icountp=0
           icountn=0
           do j=y_sta,y_end-1
              if( (sf(i,j,k,ptr_analyze+1)-zst)*(sf(i,j+1,k,20)-zst)<=0d0  )then
                 yst = yg(j) + (yg(j+1)-yg(j))/(sf(i,j+1,k,20)-sf(i,j,k,ptr_analyze+1))*(sf(i,j+1,k,20)-zst)
                 dumd= (yst-yg(j))/(yg(j+1)-yg(j))
                 chist=sf(i,j,k,ptr_analyze+2)+dumd*(sf(i,j+1,k,21)-sf(i,j,k,ptr_analyze+2))
                 ust=sf(i,j,k, 1)+dumd*(sf(i,j+1,k, 1)-sf(i,j,k, 1))
                 vst=sf(i,j,k, 2)+dumd*(sf(i,j+1,k, 2)-sf(i,j,k, 2))
                 tst=sf(i,j,k, 6)+dumd*(sf(i,j+1,k, 6)-sf(i,j,k, 6))
                 tau_mix=zst*zst/chist
                 
                 if(yst>=0d0)then
                    icountp=icountp+1
                    if(icountp==1.and.dabs(ust)>1d-5)then
                       tau_up=tau_up+1.d0/ust*dxu(i)
                    endif
                    write(333,'(8e16.7)') xg(i), yst, chist, ust, vst, tst, tau_up, tau_mix
                 else
                    icountn=icountn+1
                    if(icountn==1.and.dabs(ust)>1d-5)then
                       tau_un=tau_un+1.d0/ust*dxu(i)
                    endif
                    write(444,'(8e16.7)') xg(i), yst, chist, ust, vst, tst, tau_un, tau_mix
                 endif
              endif
           enddo
        enddo
        
     else
        write(*,*) "****  Warning *****"
        write(*,*) "      Zst profile for 3D domain is under construction "
     endif
     close(333)
     close(444)
     if(.not.flag_pv) cycle
     
     if(flag_particle)then
        write(filenumber,'(i8.8)')step
        open(ifn,file='./particle/'//'p'//filenumber//'.dat',form='unformatted')!,convert='big_endian')
        read(ifn) p_on_all
        allocate (pval_all(nv+noc,p_on_all))
        if(p_on_all>=1)then
           read(ifn) pval_all
        end if
        close(ifn)
     end if

     write(*,*)"now writing vts file"
     if(.true.)then
        call pv_output_initial( x_sta,x_end,y_sta,y_end,z_sta,z_end,ns-3,.True.)
        ! pv_output_initial( x_sta,x_end,y_sta,y_end,z_sta,z_end, number of outputed scalars, output vector or not)

        ! ptr_states
        call pv_input_scalar(sf(x_sta:x_end,y_sta:y_end,z_sta:z_end,ptr_states+1),'Density')
        call pv_input_scalar(sf(x_sta:x_end,y_sta:y_end,z_sta:z_end,ptr_states+2),'Pressure')
        call pv_input_scalar(sf(x_sta:x_end,y_sta:y_end,z_sta:z_end,ptr_states+3),'Temperature')
        call pv_input_scalar(sf(x_sta:x_end,y_sta:y_end,z_sta:z_end,ptr_states+4),'Enthaply')

         ! ptr_Y
         if (flag_Y) then
            do kk = 1, nf
               call pv_input_scalar(sf(x_sta:x_end,y_sta:y_end,z_sta:z_end,ptr_Y+kk),'Y_'//species_names_fk3(kk))
            end do
         end if

        ! ptr_analyze
        call pv_input_scalar(sf(x_sta:x_end,y_sta:y_end,z_sta:z_end,ptr_analyze+1),'Z')
        call pv_input_scalar(sf(x_sta:x_end,y_sta:y_end,z_sta:z_end,ptr_analyze+2),'chi')
        call pv_input_scalar(sf(x_sta:x_end,y_sta:y_end,z_sta:z_end,ptr_analyze+3),'Z_H')
        call pv_input_scalar(sf(x_sta:x_end,y_sta:y_end,z_sta:z_end,ptr_analyze+4),'Z_N')
        call pv_input_scalar(sf(x_sta:x_end,y_sta:y_end,z_sta:z_end,ptr_analyze+5),'H_N_ratio')
        call pv_input_scalar(sf(x_sta:x_end,y_sta:y_end,z_sta:z_end,ptr_analyze+6),'Dh')
        call pv_input_scalar(sf(x_sta:x_end,y_sta:y_end,z_sta:z_end,ptr_analyze+7),'SumY')
        call pv_input_scalar(sf(x_sta:x_end,y_sta:y_end,z_sta:z_end,ptr_analyze+8),'FI')
        call pv_input_scalar(sf(x_sta:x_end,y_sta:y_end,z_sta:z_end,ptr_analyze+9),'lambda_e')
        call pv_input_scalar(sf(x_sta:x_end,y_sta:y_end,z_sta:z_end,ptr_analyze+10),'index_maxEI')
        call pv_input_scalar(sf(x_sta:x_end,y_sta:y_end,z_sta:z_end,ptr_analyze+11),'index_maxPI')
        call pv_input_scalar(sf(x_sta:x_end,y_sta:y_end,z_sta:z_end,ptr_analyze+12),'rop_ith')
         
        ! ptr_EI
        if (flag_EI) then
            call pv_input_scalar(sf(x_sta:x_end,y_sta:y_end,z_sta:z_end,ptr_EI+1),'EI_T')
            do kk = 2, nf
               call pv_input_scalar(sf(x_sta:x_end,y_sta:y_end,z_sta:z_end,ptr_EI+kk),'EI_'//species_names_cema(kk-1))
            end do
         end if
        
        call pv_input_vector(sf(x_sta:x_end,y_sta:y_end,z_sta:z_end,1), &
                             sf(x_sta:x_end,y_sta:y_end,z_sta:z_end,2), &
                            sf(x_sta:x_end,y_sta:y_end,z_sta:z_end,3))

        call pv_input_grid(xg(x_sta:x_end),yg(y_sta:y_end),zg(z_sta:z_end))
        if(z_sta==z_end)then
           write(gridnumber,'(i5.5)')z_sta
           call pv_output_finalize(step,'./xy_plane/','pv_xy_z'//gridnumber)            
        else if(x_sta==x_end)then
           write(gridnumber,'(i5.5)')x_sta
           call pv_output_finalize(step,'./yz_plane/','pv_yz_x'//gridnumber)
        else if(y_sta==y_end)then
           write(gridnumber,'(i5.5)')y_sta
           call pv_output_finalize(step,'./zx_plane/','pv_zx_y'//gridnumber)
        else 
           call pv_output_finalize(step,'./3D/','pv3D')
        endif
        
     endif

     if(flag_particle)call pv_particle_output
     if(flag_particle)deallocate(pval_all)     
     call deallocation_cema()
     
  end do mainloop

  stop
end program pv