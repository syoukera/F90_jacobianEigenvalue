program pv
  use globals
  use pv3D
  use cema, only: calc_cema, allocation_cema, deallocation_cema
  implicit none
  double precision::dYfdx,dYfdy,dYodx,dYody

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
    read(200,'()') !--- output information ---!
    read(200,*) nf
    read(200,*) ns
    read(200,'()') !--- flags ---!
    read(200,*) flag_particle
    read(200,*) flag_Z
    read(200,*) flag_pv    
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
    write(*,*) "nf, ns  ", nf, ns
    write(*,*) "particle",flag_particle
    write(*,*) "plane   ",flag_plane

    
  call allocation
  call allocation_cema()

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

           do k=max(ksta-kbd,z_sta),min(kend+kbd,z_end)
           do j=max(jsta-jbd,y_sta),min(jend+jbd,y_end)
           do i=max(ista-ibd,x_sta),min(iend+ibd,x_end) 
              ! do k=max(ksta,z_sta),min(kend,z_end)
              ! do j=max(jsta,y_sta),min(jend,y_end)
              ! do i=max(ista,x_sta),min(iend,x_end) 
              
              sf(i,j,k,1)=u_local(i,j,k)
              sf(i,j,k,2)=v_local(i,j,k)
              sf(i,j,k,3)=w_local(i,j,k)
              sf(i,j,k,4)=r_local(i,j,k)
              sf(i,j,k,5)=p_local(i,j,k)
              sf(i,j,k,6)=t_local(i,j,k)
              sf(i,j,k,7)=h_local(i,j,k)
              
              sf(i,j,k,8)=y_local(i,j,k,1)      !N2
              sf(i,j,k,9)=y_local(i,j,k,5)      !O2
              sf(i,j,k,10)=y_local(i,j,k,7)     !OH
              sf(i,j,k,11)=y_local(i,j,k,8)     !H2
              sf(i,j,k,12)=y_local(i,j,k,9)     !H2O
              sf(i,j,k,13)=y_local(i,j,k,12)     !NH3
              sf(i,j,k,14)=y_local(i,j,k,13)     !NH2
              sf(i,j,k,15)=y_local(i,j,k,13)     !NH
              sf(i,j,k,16)=y_local(i,j,k,14)     !NO
              sf(i,j,k,17)=y_local(i,j,k,22)     !NO2
              sf(i,j,k,18)=y_local(i,j,k,23)     !N2O
              ! sf(i,j,k,19)=y_local(i,j,k,33)     !OH*

              sf(i,j,k,26)=sum(y_local(i,j,k,:))     !sumY

              dyodx=(y_local(i+1,j,k, 5)-y_local(i-1,j,k, 5))/(dxg(i)*dxg(i+1))
              dyody=(y_local(i,j+1,k, 5)-y_local(i,j-1,k, 5))/(dyg(j)*dyg(j+1))
              dyfdx=(y_local(i+1,j,k,12)-y_local(i-1,j,k,12))/(dxg(i)*dxg(i+1))
              dyfdy=(y_local(i,j+1,k,12)-y_local(i,j-1,k,12))/(dyg(j)*dyg(j+1))
              sf(i,j,k,27)=(dyodx*dyfdx+dyody*dyfdy+1d-12)/(sqrt(dyodx*dyodx+dyody*dyody)*sqrt(dyfdx*dyfdx+dyfdy*dyfdy)+1d-12)
              
              call calc_z(y_local(i,j,k,:),sf(i,j,k,20),sf(i,j,k,22),sf(i,j,k,23),sf(i,j,k,24))
              call calc_Dh(y_local(i,j,k,:),r_local(i,j,k),t_local(i,j,k),sf(i,j,k,25))
              call calc_cema(y_local(i,j,k,:),t_local(i,j,k),sf(i,j,k,28),sf(i,j,k,29))
              
           end do
           end do
           end do
           
           call deallocation
         !   call deallocation_cema
        end if
        
     end do

     sf(:,:,:,21)=0d0
     if(nz==1)then
        k=z_sta
        do j=y_sta+1,y_end-1
           do i=x_sta+1,x_end-1
              dumdx=(dxg(i)*dxg(i)*sf(i+1,j,k,20)+(dxg(i+1)*dxg(i+1)-dxg(i)*dxg(i))*sf(i,j,k,20)-dxg(i+1)*dxg(i+1)*sf(i-1,j,k,20)) &
                   / (dxg(i)*dxg(i+1)*(dxg(i)+dxg(i+1))) 
              dumdy=(dyg(j)*dyg(j)*sf(i,j+1,k,20)+(dyg(j+1)*dyg(j+1)-dyg(j)*dyg(j))*sf(i,j,k,20)-dyg(j+1)*dyg(j+1)*sf(i,j-1,k,20)) &
                   / (dyg(j)*dyg(j+1)*(dyg(j)+dyg(j+1))) 
              sf(i,j,k,21)=2d0 * sf(i,j,k,25) * (dumdx*dumdx+dumdy*dumdy)
           enddo
        enddo
     else
        do k=z_sta+1,z_end-1
           do j=y_sta+1,y_end-1
              do i=x_sta+1,x_end-1
                 dumdx=(dxg(i)*dxg(i)*sf(i+1,j,k,20)+(dxg(i+1)*dxg(i+1)-dxg(i)*dxg(i))*sf(i,j,k,20) &
                      - dxg(i+1)*dxg(i+1)*sf(i-1,j,k,20)) &
                      / (dxg(i)*dxg(i+1)*(dxg(i)+dxg(i+1))) 
                 dumdy=(dyg(j)*dyg(j)*sf(i,j+1,k,20)+(dyg(j+1)*dyg(j+1)-dyg(j)*dyg(j))*sf(i,j,k,20) &
                      - dyg(j+1)*dyg(j+1)*sf(i,j-1,k,20)) &
                      / (dyg(j)*dyg(j+1)*(dyg(j)+dyg(j+1))) 
                 dumdz=(dzg(k)*dzg(k)*sf(i,j,k+1,20)+(dzg(k+1)*dzg(k+1)-dzg(k)*dzg(k))*sf(i,j,k,20) &
                      - dzg(k+1)*dzg(k+1)*sf(i,j,k-1,20)) &
                      / (dzg(k)*dzg(k+1)*(dzg(k)+dzg(k+1))) 
                 sf(i,j,k,21)=2d0 * sf(i,j,k,25) * (dumdx*dumdx+dumdy*dumdy+dumdz*dumdz)
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
              if( (sf(i,j,k,20)-zst)*(sf(i,j+1,k,20)-zst)<=0d0  )then
                 yst = yg(j) + (yg(j+1)-yg(j))/(sf(i,j+1,k,20)-sf(i,j,k,20))*(sf(i,j+1,k,20)-zst)
                 dumd= (yst-yg(j))/(yg(j+1)-yg(j))
                 chist=sf(i,j,k,21)+dumd*(sf(i,j+1,k,21)-sf(i,j,k,21))
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
        call pv_output_initial( x_sta,x_end,y_sta,y_end,z_sta,z_end,26,.True.)
        ! pv_output_initial( x_sta,x_end,y_sta,y_end,z_sta,z_end, number of outputed scalars, output vector or not)
        call pv_input_scalar(sf(x_sta:x_end,y_sta:y_end,z_sta:z_end,4),'Density')
        call pv_input_scalar(sf(x_sta:x_end,y_sta:y_end,z_sta:z_end,5),'Pressure')
        call pv_input_scalar(sf(x_sta:x_end,y_sta:y_end,z_sta:z_end,6),'Temperature')
        call pv_input_scalar(sf(x_sta:x_end,y_sta:y_end,z_sta:z_end,7),'Enthaply')

        call pv_input_scalar(sf(x_sta:x_end,y_sta:y_end,z_sta:z_end,8),'Y_N2')
        call pv_input_scalar(sf(x_sta:x_end,y_sta:y_end,z_sta:z_end,9),'Y_O2')
        call pv_input_scalar(sf(x_sta:x_end,y_sta:y_end,z_sta:z_end,10),'Y_OH')
        call pv_input_scalar(sf(x_sta:x_end,y_sta:y_end,z_sta:z_end,11),'Y_H2')
        call pv_input_scalar(sf(x_sta:x_end,y_sta:y_end,z_sta:z_end,12),'Y_H2O')
        call pv_input_scalar(sf(x_sta:x_end,y_sta:y_end,z_sta:z_end,13),'Y_NH3')
        call pv_input_scalar(sf(x_sta:x_end,y_sta:y_end,z_sta:z_end,14),'Y_NH2')
        call pv_input_scalar(sf(x_sta:x_end,y_sta:y_end,z_sta:z_end,15),'Y_NH')
        call pv_input_scalar(sf(x_sta:x_end,y_sta:y_end,z_sta:z_end,16),'Y_NO')
        call pv_input_scalar(sf(x_sta:x_end,y_sta:y_end,z_sta:z_end,17),'Y_NO2')
        call pv_input_scalar(sf(x_sta:x_end,y_sta:y_end,z_sta:z_end,18),'Y_N2O')
        call pv_input_scalar(sf(x_sta:x_end,y_sta:y_end,z_sta:z_end,19),'Y_OHr')
        call pv_input_scalar(sf(x_sta:x_end,y_sta:y_end,z_sta:z_end,20),'Z')
        call pv_input_scalar(sf(x_sta:x_end,y_sta:y_end,z_sta:z_end,21),'chi')
        call pv_input_scalar(sf(x_sta:x_end,y_sta:y_end,z_sta:z_end,22),'Z_H')
        call pv_input_scalar(sf(x_sta:x_end,y_sta:y_end,z_sta:z_end,23),'Z_N')
        call pv_input_scalar(sf(x_sta:x_end,y_sta:y_end,z_sta:z_end,24),'H_N_ratio')
        call pv_input_scalar(sf(x_sta:x_end,y_sta:y_end,z_sta:z_end,25),'Dh')
        call pv_input_scalar(sf(x_sta:x_end,y_sta:y_end,z_sta:z_end,26),'SumY')
        call pv_input_scalar(sf(x_sta:x_end,y_sta:y_end,z_sta:z_end,27),'FI')
        call pv_input_scalar(sf(x_sta:x_end,y_sta:y_end,z_sta:z_end,28),'lambda_e')
        call pv_input_scalar(sf(x_sta:x_end,y_sta:y_end,z_sta:z_end,29),'index_maxEI')

        
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