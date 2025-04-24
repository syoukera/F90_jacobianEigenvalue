module globals
  integer(4) :: nx,ny,nz,nf,ns,nq,nrf,nrb,nrp,ista,jsta,ksta,iend,jend,kend,ibd,jbd,kbd,x_sta,x_end,y_sta,y_end,z_sta,z_end
  integer(4) :: iprocs,jprocs,kprocs,myrank,myrank_i,myrank_j,myrank_k
  integer(4) :: i,j,k,ii,jj,kk,step,ifn,is,m,n,step0,step1,step2,l,step0_avg,step1_avg,step2_avg
  double precision :: time,time0,time1,time2,ref, ugconst
  double precision , allocatable :: xu(:),xg(:),dxu(:),dxg(:)
  double precision , allocatable :: yv(:),yg(:),dyv(:),dyg(:)
  double precision , allocatable :: zw(:),zg(:),dzw(:),dzg(:)
  double precision , allocatable :: u_local(:,:,:),v_local(:,:,:),w_local(:,:,:)
  double precision , allocatable :: r_local(:,:,:),p_local(:,:,:),h_local(:,:,:),h_chem_local(:,:,:),equiv_ratio_local(:,:,:)
  double precision , allocatable :: heat_rel_local(:,:,:),mass_evap_local(:,:,:)
  double precision , allocatable :: t_local(:,:,:),c_local(:,:,:),z_local(:,:,:),rate_local(:,:,:,:)
  double precision , allocatable :: g_local(:,:,:),fmsoot_local(:,:,:),fnsoot_local(:,:,:)
  double precision , allocatable :: dp_local(:,:,:),su_local(:,:,:)
  double precision , allocatable :: rate1st2D(:,:,:),rate1st3D(:,:,:),rate2nd2D(:,:,:),rate2nd3D(:,:,:)
  double precision , allocatable :: flame_thickness(:,:,:),flame_efficiency(:,:,:)
  double precision , allocatable :: evap2D(:,:,:),evap3D(:,:,:),heat2D(:,:,:),heat3D(:,:,:)
  double precision , allocatable :: u(:,:,:),v(:,:,:),w(:,:,:),r(:,:,:),p(:,:,:),h(:,:,:),h_chem(:,:,:)
  double precision , allocatable :: temp(:,:,:),c(:,:,:),z(:,:,:),g(:,:,:),fmsoot(:,:,:),fnsoot(:,:,:)
  double precision , allocatable :: p_avg(:,:,:)
  double precision , allocatable :: pdf(:,:,:),z_val(:,:,:)
  double precision , allocatable :: y(:,:,:,:),y_local(:,:,:,:)
  double precision , allocatable :: sq_all(:,:,:,:),sf(:,:,:,:),q(:,:,:,:),sf_phase(:,:,:,:),sf_avg(:,:,:,:)
  double precision , allocatable :: ccpl(:,:), ccph(:,:), tmpr(:,:)
  double precision , allocatable , save :: cmu(:,:),clam(:,:)
  double precision , allocatable , save :: cdif(:,:,:)
  character :: filenumber*8,cpunumber*5,filenumber0*8,filenumber1*8,filenumber2*8,reaction_number*2, gridnumber*5

  integer(4) , parameter :: p_total =1000000,split=5,noc=1,nv=14+1
  integer(4) ,      save :: p_on,substep,npaa,p_on_all
  !integer(4) ,                      save :: p_in,igx,igy,igz
  logical, dimension(p_total) :: p_label=.false.
  double precision                     :: mass_in,mass_rate,equiv_ratio,wg_avg,p_dmin,interval=0.d0,heat2Dall,heat3Dall
  double precision, dimension(p_total) :: p_x,p_u,gp_u0
  double precision, dimension(p_total) :: p_y,p_v,gp_v0
  double precision, dimension(p_total) :: p_z,p_w,gp_w0
  double precision, dimension(p_total) :: p_r,p_mass,p_d,p_ds ,p_temp,cl,p_h
  double precision, dimension(p_total) :: gp_mass0,gp_temp0
  double precision, dimension(p_total) :: p_cl,p_id,parcel_on
  double precision, dimension(p_total,noc) , save :: c_mass

  double precision :: integral_time,previous_time

  double precision, dimension(nv+noc,p_total) :: pval
  double precision, allocatable , save :: pval_all(:,:), pval_all_tmp(:,:)

  integer(4) :: rows, step_list0, step_list1, step_list2
  integer,allocatable :: step_list(:)
  double precision,allocatable :: time_list(:)

  ! double precision , dimension(1-ibd:nx+4,1-ibd:ny+4,1-ibd:nz+4,1:ns)::sf

  double precision , parameter :: na=6.02d+26

  integer(4) , dimension(1:100,1:8) , save :: step_phase=0

  double precision :: phi
!  double precision , dimension(1:6) , save :: wg
  double precision , allocatable :: wg(:)

  double precision , allocatable , save       :: ri(:,:,:),rayleighindex(:,:,:)

  logical, save :: flag_particle
  logical, save :: flag_Z,flag_pv,flag_Y,flag_EI
  double precision , allocatable , save       :: coeffZ(:,:)
  double precision :: zh1,zh2,zc1,zc2,zn1,zn2,zo1,zo2,owc,owh,own,owo,cfuel,hfuel,nfuel,ofuel,nu_o2,zc,zh,zn,zo, wc,wh,wo,wn
  double precision :: dumd,dumdx,dumdy,dumdz,xo2,zst,chist,Tst,ust,vst,yst,tau_mix,tau_up,tau_un
  integer::icountp,icountn
  
  integer,save :: flag_plane 
!================production & destruction=========!ayukawa
  double precision :: h_reac_in,e_reac_in,r_reac_in,p_reac_in,temp0,dum_dble,dum_dble2
  double precision :: dt=1.d-7
  real(4), allocatable :: array3d(:,:,:,:)
  
  character(99) , save :: chem_dir
  character species_names_fk3(1000)*24

end module globals

module pv3D
  !===============================
  !    2017/10/16 Updated
  !===============================
  !     MEMO
  ! paraview可視化ファイル作成用モジュール
  ! binary形式での書き出し。読み込みが早い。
  ! float64にしているが、float32にすれば
  ! 容量を半減させることができるので、それも
  ! 選択式に出来たらよい。
  !
  !     Data format of AppendedData
  ! _N1[Data1]N2[Data2]...
  ! 冒頭の"_"は必須。N1,N2はデータのバイト数。
  ! Data1の総バイト数を記述。offsetは、N1,N2の冒頭部分を指定。
  !
  ! Time Stamp Output追加 !17/10/16
  !
  !
  !
  !===============================


  !===============================
  !     Subroutine's Scopes
  !===============================
  ! public  :: pv_output, region_input

  !===============================
  !       Structure Type
  !===============================
  ! type vector
  !   double precision :: x,y,z
  ! end type vector
  ! type scalar_value
  !   character(len=20)     :: name
  !   double precision , allocatable :: array(:,:,:)
  ! end type scalar_value

  !===============================
  !     Parameter Settings
  !===============================
  integer(4) , private , parameter :: ifn=10
  integer(4) , private :: ifnwrite
  integer(4) , private , parameter :: bytes=4 ! 8 or 4
  integer(4) , private , save      :: num_output,num_output_pointer
  double precision , private , save , allocatable      :: output_scalar_all(:,:,:,:)
  real , private , save , allocatable      :: output_scalar_all4(:,:,:,:) !yada
  logical , private , save      :: vecflag
  integer(4) , private , parameter :: ns=11
  integer(4) , private , parameter :: nf=9
  ! integer(4) , private , parameter :: step0
  ! integer(4) , private , parameter :: step1
  ! integer(4) , private , parameter :: step2
  ! integer(4) , private , parameter :: nx
  ! integer(4) , private , parameter :: ny
  ! integer(4) , private , parameter :: nz
  ! integer(4) , private , parameter :: ibd
  ! integer(4) , private , parameter :: jbd
  ! integer(4) , private , parameter :: kbd
  integer(4) , private :: nq,ista,jsta,ksta,iend,jend,kend
  integer(4) , private , save :: iprocs,jprocs,kprocs,myrank,myrank_i,myrank_j,myrank_k
  character(len=8), private :: data_type


  !===============================
  !          Variables
  !===============================
  integer(4) , private :: i,j,k,ii,jj,kk
  integer(4) , private :: x_sta_reg,x_end_reg,y_sta_reg,y_end_reg,z_sta_reg,z_end_reg
  integer(4) , private :: grid_offset, vector_offset
  integer(4) , private :: data_size
  integer(4) , private :: block_data_size
  integer(4) , private , allocatable :: scalar_offset(:)
  character(8) , private :: filenumber, gridnumber
  character(16),allocatable , private , save :: str_varname(:)
  character(16),allocatable , private :: str_format(:)
  double precision , private , allocatable , save :: u_local(:,:,:),v_local(:,:,:),w_local(:,:,:)
  double precision , private , allocatable :: r_local(:,:,:)
  double precision , private , allocatable :: p_local(:,:,:)
  double precision , private , allocatable :: t_local(:,:,:)
  double precision , private , allocatable :: h_local(:,:,:)
  double precision , private , allocatable :: y_local(:,:,:,:)
  double precision , private , allocatable :: z_local(:,:,:)
  double precision , private , allocatable :: c_local(:,:,:)
  double precision , private , allocatable :: g_local(:,:,:)
  double precision , private , allocatable :: fmsoot_local(:,:,:),fnsoot_local(:,:,:)
  double precision , private :: time
  double precision , private , allocatable , save :: xg(:)!,dxu(:),dxg(:)
  double precision , private , allocatable , save :: yg(:)!,dyv(:),dyg(:)
  double precision , private , allocatable , save :: zg(:)!,dzw(:),dzg(:
  real(4) , private , allocatable:: array3d(:,:,:,:)

  ! type(scalar_value) , allocatable , private  :: valuebox(:)

contains

  !============================!
  !=== DATA INPUT FUNCTIONS ===!
  !============================!
  subroutine pv_output_initial(x_sta_in,x_end_in,y_sta_in,y_end_in,z_sta_in,z_end_in,scalar_num,vector_flag)
    implicit none
    integer(4) , intent(in) :: x_sta_in,x_end_in,y_sta_in,y_end_in,z_sta_in,z_end_in
    integer(4) , intent(in) :: scalar_num
    logical , intent(in) :: vector_flag

    integer(4) :: l

    num_output=scalar_num
    vecflag=vector_flag

    select case(bytes)
    case(4)
      data_type='Float32'
    case(8)
      data_type='Float64'
    end select

    x_sta_reg=x_sta_in
    x_end_reg=x_end_in
    y_sta_reg=y_sta_in
    y_end_reg=y_end_in
    z_sta_reg=z_sta_in
    z_end_reg=z_end_in

    allocate(output_scalar_all(x_sta_reg:x_end_reg,y_sta_reg:y_end_reg,z_sta_reg:z_end_reg,1:num_output))
    allocate(output_scalar_all4(x_sta_reg:x_end_reg,y_sta_reg:y_end_reg,z_sta_reg:z_end_reg,1:num_output))
    allocate(u_local(x_sta_reg:x_end_reg,y_sta_reg:y_end_reg,z_sta_reg:z_end_reg))
    allocate(v_local(x_sta_reg:x_end_reg,y_sta_reg:y_end_reg,z_sta_reg:z_end_reg))
    allocate(w_local(x_sta_reg:x_end_reg,y_sta_reg:y_end_reg,z_sta_reg:z_end_reg))
    allocate(str_varname(1:num_output))
    allocate(scalar_offset(1:num_output))
    allocate(xg(x_sta_reg:x_end_reg),yg(y_sta_reg:y_end_reg),zg(z_sta_reg:z_end_reg))
    allocate( array3d(1:3,x_sta_reg:x_end_reg,y_sta_reg:y_end_reg,z_sta_reg:z_end_reg))

    !=== offset setting ===!
    block_data_size = size(xg)*size(yg)*size(zg)*bytes

    grid_offset=0

    scalar_offset=0
    scalar_offset(1)=grid_offset + block_data_size*3+8
    if(num_output>1)then
      do l=2,num_output
        scalar_offset(l)=scalar_offset(l-1)+block_data_size+8
      enddo
    endif

    vector_offset=0
    if(vecflag)then
      vector_offset=scalar_offset(num_output)+block_data_size+8
      data_size=vector_offset+block_data_size*3
    else
      data_size=scalar_offset(num_output)+block_data_size
    endif
    ! write(*,*) data_size , "DATA SIZE !!!"

    num_output_pointer=1

  end subroutine pv_output_initial

  subroutine pv_input_scalar(input,scalarname)
    implicit none
    double precision , intent(in)  :: input( x_sta_reg:x_end_reg,y_sta_reg:y_end_reg,z_sta_reg:z_end_reg)
    character(*) , intent(in) :: scalarname

    output_scalar_all(:,:,:,num_output_pointer)=input(:,:,:)
    str_varname(num_output_pointer)=trim(adjustl(scalarname))

    num_output_pointer=num_output_pointer+1

  end subroutine pv_input_scalar

  subroutine pv_input_vector(u,v,w)
    implicit none
    double precision , intent(in)  :: u( x_sta_reg:x_end_reg,y_sta_reg:y_end_reg,z_sta_reg:z_end_reg)
    double precision , intent(in)  :: v( x_sta_reg:x_end_reg,y_sta_reg:y_end_reg,z_sta_reg:z_end_reg)
    double precision , intent(in)  :: w( x_sta_reg:x_end_reg,y_sta_reg:y_end_reg,z_sta_reg:z_end_reg)

    u_local=u
    v_local=v
    w_local=w

  end subroutine pv_input_vector

  subroutine pv_input_grid(xgrid,ygrid,zgrid)
    implicit none
    double precision , intent(in)  :: xgrid(x_sta_reg:x_end_reg),ygrid(y_sta_reg:y_end_reg),zgrid(z_sta_reg:z_end_reg)

    xg=xgrid
    yg=ygrid
    zg=zgrid

  end subroutine pv_input_grid

  !==============================!
  !=== DATA WRITING FUNCTIONS ===!
  !==============================!
  subroutine header_output
    implicit none

    write(ifn) '<VTKFile type="StructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">'
    write(ifn) &
    '<StructuredGrid WholeExtent="',&
    trim(adjustl(int_to_char(x_sta_reg)))," ",&
    trim(adjustl(int_to_char(x_end_reg)))," ",&
    trim(adjustl(int_to_char(y_sta_reg)))," ",&
    trim(adjustl(int_to_char(y_end_reg)))," ",&
    trim(adjustl(int_to_char(z_sta_reg)))," ",&
    trim(adjustl(int_to_char(z_end_reg))),'">'
    write(ifn) &
    '<Piece Extent="',&
    trim(adjustl(int_to_char(x_sta_reg)))," ",&
    trim(adjustl(int_to_char(x_end_reg)))," ",&
    trim(adjustl(int_to_char(y_sta_reg)))," ",&
    trim(adjustl(int_to_char(y_end_reg)))," ",&
    trim(adjustl(int_to_char(z_sta_reg)))," ",&
    trim(adjustl(int_to_char(z_end_reg))),'">'

  end subroutine header_output

  subroutine footer_output
    implicit none

    write(ifn) '</Piece>'
    write(ifn) '</StructuredGrid>'
    
    write(*,*) "before call 32", bytes
    write(*,*) "x_sta_reg",x_sta_reg 
    write(*,*) "x_end_reg",x_end_reg 
    write(*,*) "y_sta_reg",y_sta_reg 
    write(*,*) "y_end_reg",y_end_reg 
    write(*,*) "z_sta_reg",z_sta_reg 
    write(*,*) "z_end_reg",z_end_reg 
 
    select case(bytes)
    case(4)
      call appenddata_output_float32
    case(8)
      call appenddata_output_float64
    end select

    write(ifn) '</VTKFile>'

  end subroutine footer_output

  subroutine pointdata_header_output
    implicit none

    write(ifn) '<PointData Scalars="'//str_varname(1)//'" Vectors="Velocity">'

  end subroutine pointdata_header_output

  subroutine pointdata_footer_output
    implicit none

    write(ifn) '</PointData>'

  end subroutine pointdata_footer_output

  subroutine scalarall_output
    implicit none

    integer(4) :: l

    do l=1,num_output
      call scalar_output(output_scalar_all(:,:,:,l),str_varname(l),l)
    enddo

  end subroutine scalarall_output

  subroutine scalar_output(output_array,scalarname,l)
    implicit none

    double precision :: output_array(x_sta_reg:x_end_reg,y_sta_reg:y_end_reg,z_sta_reg:z_end_reg)
    integer(4) , intent(in) :: l
    character(*) , intent(in) :: scalarname


    write(ifn) &
    '<DataArray type="'//trim(adjustl(data_type))//&
    '" Name="'// trim(adjustl(scalarname)) //'" RangeMin="',&
    trim(adjustl(flt_to_char(minval(output_array)))),&
    '" RangeMax="',&
    trim(adjustl(flt_to_char(maxval(output_array)))),&
    '" format="appended" offset="',&
    trim(adjustl(int_to_char8(scalar_offset(l)))),&
    '" />'

  end subroutine scalar_output

  subroutine vector_output
    implicit none

    write(ifn) &
    '<DataArray type="'//trim(adjustl(data_type))//&
    '" Name="Velocity" NumberOfComponents="3" RangeMin="',&
    ! trim(adjustl(flt_to_char(min(minval(u_local),minval(v_local),minval(w_local))))),&
    '" RangeMax="',&
    ! trim(adjustl(flt_to_char(max(maxval(u_local),maxval(v_local),maxval(w_local))))),&
    '" format="appended" offset="',&
    trim(adjustl(int_to_char8(vector_offset))),&
    '" />'

  end subroutine vector_output

  subroutine griddata_output
    implicit none

    write(ifn) '<Points>'
    write(ifn) &
    '<DataArray type="'//trim(adjustl(data_type))//&
    '" Name="Points" NumberOfComponents="3" format="appended" RangeMin="" RangeMax="" offset="',&
    trim(adjustl(int_to_char8(grid_offset))),&
    '" />'
    write(ifn) '</Points>'

  end subroutine griddata_output

  subroutine appenddata_output_float64
    implicit none

    integer :: l
    ! integer :: output_data_size=0
    double precision :: array3d(1:3,x_sta_reg:x_end_reg,y_sta_reg:y_end_reg,z_sta_reg:z_end_reg)

    do k=z_sta_reg,z_end_reg
      do j=y_sta_reg,y_end_reg
        do i=x_sta_reg,x_end_reg
          array3d(1,i,j,k)=xg(i)
          array3d(2,i,j,k)=yg(j)
          array3d(3,i,j,k)=zg(k)
        enddo
      enddo
    enddo

    write(ifn) '<AppendedData encoding="raw">'
    write(ifn) '_'
    write(ifn) data_size


    !=== Grid data output ===!
    ! do k=z_sta_reg,z_end_reg
    ! do j=y_sta_reg,y_end_reg
    ! do i=x_sta_reg,x_end_reg
    ! write(ifn) xg(i),yg(j),zg(k)
    write(ifn) array3d
    ! write(*,'(3f20.12)') xg(i),yg(j),zg(k)
    ! enddo
    ! enddo
    ! enddo
    ! output_data_size=output_data_size+&
    ! size(xg(x_sta_reg:x_end_reg))*&
    ! size(yg(y_sta_reg:y_end_reg))*&
    ! size(zg(z_sta_reg:z_end_reg))*8

    !=== Scalar data output ===!
    do l=1,num_output
      write(ifn) block_data_size
      ! do k=z_sta_reg,z_end_reg
      ! do j=y_sta_reg,y_end_reg
      ! do i=x_sta_reg,x_end_reg
      ! if(isnan(output_scalar_all(i,j,k,l)))then
      ! write(ifn) dble(0.d0)
      ! else
      ! if(dabs(output_scalar_all(i,j,k,l))<=1.d-30)then
      ! write(ifn) dble(0.d0)
      ! else
      write(ifn) output_scalar_all(:,:,:,l)
      ! endif
      ! endif
      ! enddo
      ! enddo
      ! enddo

    enddo

    !=== Vector data output ===!
    if(vecflag)then
      do k=z_sta_reg,z_end_reg
        do j=y_sta_reg,y_end_reg
          do i=x_sta_reg,x_end_reg
            array3d(1,i,j,k)=u_local(i,j,k)
            array3d(2,i,j,k)=v_local(i,j,k)
            array3d(3,i,j,k)=w_local(i,j,k)
          enddo
        enddo
      enddo

      write(ifn) block_data_size*3
      write(ifn) array3d

    endif

    write(ifn) '</AppendedData>'


  end subroutine appenddata_output_float64

  subroutine appenddata_output_float32
    implicit none

    integer :: l
    ! integer :: output_data_size=0
!    real(4) :: array3d(1:3,x_sta_reg:x_end_reg,y_sta_reg:y_end_reg,z_sta_reg:z_end_reg)

    write(*,*) "************ 1 ************* "

    do k=z_sta_reg,z_end_reg
      do j=y_sta_reg,y_end_reg
        do i=x_sta_reg,x_end_reg
          array3d(1,i,j,k)=real(xg(i))
          array3d(2,i,j,k)=real(yg(j))
          array3d(3,i,j,k)=real(zg(k))
          do l=1,num_output
             output_scalar_all4(i,j,k,l)=real(output_scalar_all(i,j,k,l))
          end do
       enddo
      enddo
    enddo

    write(*,*)"array3d.size=",size(array3d)
    write(*,*) "************ 2 ************* "

    write(ifn) '<AppendedData encoding="raw">'
    write(ifn) '_'
    write(ifn) data_size
    write(*,*)"datasize=",data_size

    write(*,*) "************ 3 ************* "

    !=== Grid data output ===!
    write(ifn) array3d

    write(*,*) "************ 4 ************* "

    !=== Scalar data output ===!
    write(*,*)"num_output",num_output
    do l=1,num_output
      write(ifn) block_data_size
      write(*,*)"outputscalar_all.size=",size(output_scalar_all4)
      write(*,*)"block_data_size",block_data_size
      !write(ifn) real(output_scalar_all(:,:,:,l))
      write(ifn) output_scalar_all4(:,:,:,l)
      write(*,*)"finished! l=",l
    enddo

    write(*,*) "************ 5 ************* "

    !=== Vector data output ===!
    if(vecflag)then
      do k=z_sta_reg,z_end_reg
        do j=y_sta_reg,y_end_reg
          do i=x_sta_reg,x_end_reg
            array3d(1,i,j,k)=real(u_local(i,j,k))
            array3d(2,i,j,k)=real(v_local(i,j,k))
            array3d(3,i,j,k)=real(w_local(i,j,k))
          enddo
        enddo
      enddo

      write(ifn) block_data_size*3
      write(ifn) array3d

    endif

    write(*,*) "************ 6 ************* "

    write(ifn) '</AppendedData>'


  end subroutine appenddata_output_float32

  subroutine pv_output_finalize(step,dir,filename)
    implicit none
    integer(4) , intent(in) :: step
    character(*) , intent(in) :: dir,filename
    character(10) :: filenumber

    if(num_output_pointer/=num_output+1) then
      write(*,*) 'num_output not correct'
    else
      write(filenumber,'(i10.10)') step

      open(ifn,file=dir//'/'//filename//'_'//filenumber//'.vts',status='replace',access='stream',form='unformatted')

      call header_output


      call pointdata_header_output
      if(vecflag) call vector_output
      call scalarall_output
      call pointdata_footer_output
      call griddata_output
      call footer_output

      close(ifn)

    endif

    deallocate(output_scalar_all)
    deallocate(output_scalar_all4)
    deallocate(u_local)
    deallocate(v_local)
    deallocate(w_local)
    deallocate(str_varname)
    deallocate(scalar_offset)
    deallocate(xg,yg,zg)
    deallocate(array3d)

  end subroutine pv_output_finalize

  !--- time stamp output ---
  subroutine pv_time_initialize


    ifnwrite=99
    open(ifnwrite,file='timestep.pvd',status='replace')

  end subroutine pv_time_initialize

  subroutine pv_time_finalize
    implicit none

    close(ifnwrite)

  end subroutine pv_time_finalize

  subroutine pv_time_header_output
    implicit none

    write(ifnwrite,'(A)') '<?xml version="1.0"?>'
    write(ifnwrite,'(A)') '<VTKFile type="Collection" version="0.1">'
    write(ifnwrite,'(A)') '<Collection>'

  end subroutine pv_time_header_output

  subroutine pv_time_footer_output
    implicit none

    write(ifnwrite,'(A)') '</Collection>'
    write(ifnwrite,'(A)') '</VTKFile>'
    close(ifnwrite)

  end subroutine pv_time_footer_output

  subroutine timestep_output(time,stepnumber)
    implicit none

    double precision , intent(in) :: time
    character(10)  :: file_preword,file_postword
    ! character(10)  :: file_preword_particle,file_postword_particle
    character(8) , intent(in) :: stepnumber

    file_preword="pv2D_"
    file_postword=".vts"
    write(*,*) 'writing ' // file_preword // stepnumber // trim(adjustl(file_postword))
    write(ifnwrite,'(A,1e20.10,A)') '<DataSet timestep="' , time ,&
    '" group="" part="0" file="'//trim(adjustl(file_preword))//&
    trim(adjustl(stepnumber))//trim(adjustl(file_postword))//'" />'

  end subroutine timestep_output

  character(100) function int_to_char(input)
    integer(4) , intent(in) :: input

    write(int_to_char,'(1I100)') input

  end function int_to_char


  character(100) function int_to_char8(input)
    integer(4) , intent(in) :: input

    write(int_to_char8,'(1I100)') input

  end function int_to_char8


  character(100) function flt_to_char(input)
    double precision , intent(in) :: input

    write(flt_to_char,'(1e12.4)') input

  end function flt_to_char

end module pv3D

module pv_particle
  use globals
  implicit none
  ! include 'mpif.h'
  private


  public :: pv_output_initialize
  public :: pv_input_scalar
  public :: pv_input_points
  public :: pv_input_velocity
  public :: pv_output_file
  public :: pv_output_finalize


  !--- Global Variables ---
  integer(4) :: step_start, step_end, step_int
!  integer(4) :: x_sta,x_end,y_sta,y_end,z_sta,z_end
  integer(4) :: pv_p_num=3

  !--- Flags ---
  logical :: filter_flag,pv_output_flag,pdf_output_flag,particle_read_flag,euler_read_flag

  !--- Particle's Variables ---
  logical , allocatable :: pval_flag(:)

  integer(4) , save :: num_output,num_output_pointer
  integer(4) , save :: p_on_all_filtered
  integer(4) , save :: p_sum,p_on_localmax,p_counter
  integer(4) , allocatable :: p_i(:),p_j(:),p_k(:)

  double precision , allocatable :: output_scalar_all(:,:),velocity_all(:,:),points_all(:,:)

  character(16),allocatable :: str_varname(:)
  character(16),allocatable :: str_format(:)

  !--- Eulerial Data's Variables ---
  integer(4) :: euler_data_number
  double precision :: xg_sta,yg_sta,zg_sta,xg_end,yg_end,zg_end

  double precision , allocatable :: all_data(:,:,:,:)

  double precision , allocatable , save :: p_d_output(:,:),p_d_output_local(:,:),p_d_all(:),p_d_all_all(:)

contains


  subroutine grid_setting
    implicit none

    xg_sta=xg(x_sta)
    yg_sta=yg(y_sta)
    zg_sta=zg(z_sta)
    xg_end=xg(x_end)
    yg_end=yg(y_end)
    zg_end=zg(z_end)
  end subroutine grid_setting


  subroutine grid_filter
    implicit none

    integer(4) :: logical_number
    pval_flag=.false.
    logical_number=0

    do ii=1,p_on_all

       if(points_all(1,ii)>=xg_sta .and. points_all(1,ii)<=xg_end)then
        if(points_all(2,ii)>=yg_sta .and. points_all(2,ii)<=yg_end)then
          if(points_all(3,ii)>=zg_sta .and. points_all(3,ii)<=zg_end)then
            pval_flag(ii)=.true.
            logical_number=logical_number+1
          endif
        endif
      endif

   enddo

    write(*,*) logical_number
    p_on_all_filtered=logical_number
  end subroutine grid_filter


  subroutine pv_output_initialize(numberofpoints,numberofscalaroutput)
    implicit none

    integer(4) , intent(in) :: numberofpoints,numberofscalaroutput

    num_output=numberofscalaroutput

    allocate(output_scalar_all(1:num_output,1:numberofpoints))
    allocate(str_varname(1:num_output))
    allocate(velocity_all(1:3,1:numberofpoints))
    allocate(points_all(1:3,1:numberofpoints))
    allocate(pval_flag(1:numberofpoints))
    pval_flag=.true.

    num_output_pointer=1

  end subroutine pv_output_initialize


  
  subroutine pv_output_file(outputdir,outputhead)
    implicit none

    character(*) ,intent(in) :: outputdir,outputhead

    open(ifn,file=outputdir//'/'//outputhead//filenumber//'.vtu',form='formatted')

    call header_output
    call variableall_output
    write(ifn,'(A)') '<CellData></CellData>'
    call Pointsdata_output
    write(ifn,'(A)') '<Cells>&
    <DataArray type="Int64" Name="connectivity" format="ascii" RangeMin="1e+299" RangeMax="-1e+299">&
    </DataArray>&
    <DataArray type="Int64" Name="offsets" format="ascii" RangeMin="1e+299" RangeMax="-1e+299">&
    </DataArray>&
    </Cells>'
    call footer_output

    close(ifn)

  end subroutine pv_output_file


  subroutine pv_output_finalize
    implicit none


    deallocate(output_scalar_all)
    deallocate(str_varname)
    deallocate(velocity_all)
    deallocate(points_all)
    deallocate(pval_flag)

  end subroutine pv_output_finalize

  subroutine header_output
    implicit none

    write(ifn,'(A)') '<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">'
    write(ifn,'(A)') '<UnstructuredGrid>'
    write(ifn,'(A,1(I8,x),A)') '<Piece NumberOfPoints="',p_on_all,'" NumberOfCells="0" >'

  end subroutine header_output

  subroutine pointdata_header_output
    implicit none

    write(ifn,'(A)') '<PointData Scalars="'//str_varname(1)//'" Vectors="Velocity" >'

  end subroutine pointdata_header_output


  subroutine pointdata_footer_output
    implicit none

    write(ifn,'(A)') '</PointData>'

  end subroutine pointdata_footer_output

  subroutine pv_input_velocity(p_u,p_v,p_w)
    implicit none

    double precision , intent(in) :: p_u(1:p_on_all),p_v(1:p_on_all),p_w(1:p_on_all)

    velocity_all(1,:)=p_u(:)
    velocity_all(2,:)=p_v(:)
    velocity_all(3,:)=p_w(:)

  end subroutine pv_input_velocity

  subroutine pv_input_points(p_x,p_y,p_z)
    implicit none

    double precision , intent(in) :: p_x(1:p_on_all),p_y(1:p_on_all),p_z(1:p_on_all)

    points_all(1,:)=p_x(:)
    points_all(2,:)=p_y(:)
    points_all(3,:)=p_z(:)

  end subroutine pv_input_points


  subroutine pv_input_scalar(input,scalarname)
    implicit none

    character(*) scalarname
    double precision , intent(in), dimension(1:p_on_all) :: input

    write(*,*) scalarname//' input scalar' , maxval(input),minval(input)

    output_scalar_all(num_output_pointer,:)=input(:)
    str_varname(num_output_pointer)=trim(adjustl(scalarname))

    write(*,*) str_varname(num_output_pointer)//' inputed scalar' , &
    maxval(output_scalar_all(num_output_pointer,:)),minval(output_scalar_all(num_output_pointer,:))

    num_output_pointer=num_output_pointer+1

  end subroutine pv_input_scalar


  subroutine scalar_output(output_array,scalarname)
    implicit none

    double precision , intent(in) :: output_array(1:p_on_all)
    double precision :: write_array(1:p_on_all)
    character(*) , intent(in) :: scalarname

    write_array=output_array!pack(output_array,pval_flag)

    write(ifn,'(A,1e16.8,A,1e16.8,A)') '<DataArray type="Float32" Name="' &
    // scalarname //'" format="ascii" RangeMin="',minval(write_array), &
    '" RangeMax="',maxval(write_array),'">'
    write(ifn,'(6e20.10)') write_array
    write(ifn,'(A)') '</DataArray>'

  end subroutine scalar_output

  subroutine variableall_output
    implicit none

    integer(4) :: l

    call pointdata_header_output
    write(*,*) num_output_pointer,num_output
    if(num_output_pointer/=num_output+1) stop 'error! not correct num'
    do l=1,num_output
      call scalar_output(output_scalar_all(l,:),str_varname(l))
    enddo
    call Vectordata_output
    call pointdata_footer_output

  end subroutine variableall_output


  subroutine Pointsdata_output
    implicit none

    double precision :: range_max,range_min

    range_max=max(maxval(points_all(1,:)),maxval(points_all(2,:)),maxval(points_all(3,:)))
    range_min=min(minval(points_all(1,:)),minval(points_all(2,:)),minval(points_all(3,:)))

    write(ifn,'(A)') '<Points>'
    write(ifn,'(A,1e16.8,A,1e16.8,A)') &
    '<DataArray type="Float32" Name="Points" NumberOfComponents="3" format="ascii" RangeMin="', &
    range_min,'" RangeMax="',range_max,'">'
    do l=1,p_on_all
      if(pval_flag(l)) write(ifn,'(3e20.10)') points_all(1,l),points_all(2,l),points_all(3,l)
    enddo
    write(ifn,'(A)') '</DataArray>'
    write(ifn,'(A)') '</Points>'

  end subroutine Pointsdata_output


  subroutine Vectordata_output
    implicit none

    double precision :: range_max,range_min

    range_max=max(maxval(velocity_all(1,:)),maxval(velocity_all(2,:)),maxval(velocity_all(3,:)))
    range_min=min(minval(velocity_all(1,:)),minval(velocity_all(2,:)),minval(velocity_all(3,:)))

    write(ifn,'(A,1e16.8,A,1e16.8,A)') &
    '<DataArray type="Float32" Name="Velocity" NumberOfComponents="3" format="ascii" RangeMin="',&
    range_min,'" RangeMax="',range_max,'">'
    do l=1,p_on_all
      if(pval_flag(l)) write(ifn,'(3e20.10)') velocity_all(1,l),velocity_all(2,l),velocity_all(3,l)
    enddo
    write(ifn,'(A)') '</DataArray>'

  end subroutine Vectordata_output

  subroutine footer_output
    implicit none!

    write(*,*) 'test'
    write(ifn,'(A)') '</Piece>'
    write(ifn,'(A)') '</UnstructuredGrid>'
    write(ifn,'(A)') '</VTKFile>'!

  end subroutine footer_output


  subroutine output_term
    implicit none

    write(*,*)'x',maxval(p_x,p_label),minval(p_x,p_label)
    write(*,*)'y',maxval(p_y,p_label),minval(p_y,p_label)
    write(*,*)'z',maxval(p_z,p_label),minval(p_z,p_label)
    write(*,*)'u',maxval(p_u,p_label),minval(p_u,p_label)
    write(*,*)'v',maxval(p_v,p_label),minval(p_v,p_label)
    write(*,*)'w',maxval(p_w,p_label),minval(p_w,p_label)
    write(*,*)'d',maxval(p_d,p_label),minval(p_d,p_label)
    write(*,*)'r',maxval(p_r,p_label),minval(p_r,p_label)
    write(*,*)'temp',maxval(p_temp,p_label),minval(p_temp,p_label)
    write(*,*)'mass',maxval(p_mass,p_label),minval(p_mass,p_label)
    write(*,*)'mass0',maxval(gp_mass0,p_label),minval(gp_mass0,p_label)

  end subroutine output_term

end module pv_particle

subroutine pv_particle_output
use globals
use pv_particle
implicit none

      call pv_output_initialize(p_on_all,3)
      call pv_input_scalar(pval_all(7,:),'P_Diameter')
!test
     if(.false.)then
      open(ifn,file='p_diameter.dat',form='formatted',position='append')!,convert='big_endian')
      write(ifn,*)pval_all(7,:)
      close(ifn)
     endif

      call pv_input_scalar(pval_all(9,:),'P_Temperature')
      call pv_input_scalar(pval_all(10,:),'P_Mass')
      call pv_input_points(pval_all(1,:),pval_all(2,:),pval_all(3,:))
      call pv_input_velocity(pval_all(4,:),pval_all(5,:),pval_all(6,:))
      if(flag_plane==1)then
      call pv_output_file('./xy_plane','pvp') !change
      else if(flag_plane==2)then
      call pv_output_file('./yz_plane','pvp') !change
      else if(flag_plane==3)then
      call pv_output_file('./zx_plane','pvp') !change
      else
      call pv_output_file('./3D','pvp') !change         
      endif
     call pv_output_finalize

end subroutine
  

subroutine allocation
  use globals

  allocate(xu(1-ibd:nx+ibd))
  allocate(xg(1-ibd:nx+ibd))
  allocate(dxu(1-ibd:nx+ibd))
  allocate(dxg(1-ibd:nx+ibd))
  allocate(yv(1-jbd:ny+jbd))
  allocate(yg(1-jbd:ny+jbd))
  allocate(dyv(1-jbd:ny+jbd))
  allocate(dyg(1-jbd:ny+jbd))
  allocate(zw(1-kbd:nz+kbd))
  allocate(zg(1-kbd:nz+kbd))
  allocate(dzw(1-kbd:nz+kbd))
  allocate(dzg(1-kbd:nz+kbd))
  !allocate(u(1-ibd:nx+ibd,1-jbd:ny+jbd,1-kbd:nz+kbd))
  !allocate(v(1-ibd:nx+ibd,1-jbd:ny+jbd,1-kbd:nz+kbd))
  !allocate(w(1-ibd:nx+ibd,1-jbd:ny+jbd,1-kbd:nz+kbd))
  !allocate(r(1-ibd:nx+ibd,1-jbd:ny+jbd,1-kbd:nz+kbd))
  !allocate(p(1-ibd:nx+ibd,1-jbd:ny+jbd,1-kbd:nz+kbd))
  !allocate(h(1-ibd:nx+ibd,1-jbd:ny+jbd,1-kbd:nz+kbd))
  !allocate(h_chem(1-ibd:nx+ibd,1-jbd:ny+jbd,1-kbd:nz+kbd))
  !allocate(temp(1-ibd:nx+ibd,1-jbd:ny+jbd,1-kbd:nz+kbd))
  !allocate(c(1-ibd:nx+ibd,1-jbd:ny+jbd,1-kbd:nz+kbd))
  !allocate(z(1-ibd:nx+ibd,1-jbd:ny+jbd,1-kbd:nz+kbd))
  !allocate(g(1-ibd:nx+ibd,1-jbd:ny+jbd,1-kbd:nz+kbd))
  !allocate(fmsoot(1-ibd:nx+ibd,1-jbd:ny+jbd,1-kbd:nz+kbd))
  !allocate(fnsoot(1-ibd:nx+ibd,1-jbd:ny+jbd,1-kbd:nz+kbd))

  !allocate(p_avg(1-ibd:nx+ibd,1-jbd:ny+jbd,1-kbd:nz+kbd))
  !allocate(pdf(1-ibd:nx+ibd,1-jbd:ny+jbd,1-kbd:nz+kbd))
  !allocate(z_val(1-ibd:nx+ibd,1-jbd:ny+jbd,1-kbd:nz+kbd))
  !allocate(y(1-ibd:nx+ibd,1-jbd:ny+jbd,1-kbd:nz+kbd,1:nf))
  !allocate(sq_all(1-ibd:nx+ibd,1-jbd:ny+jbd,1-kbd:nz+kbd,1:nq))
  !allocate(sf(1-ibd:nx+ibd,1-jbd:ny+jbd,1-kbd:nz+kbd,1:ns))
  allocate(sf(x_sta:x_end,y_sta:y_end,z_sta:z_end,1:ns))
  ! allocate(sf_phase(x_sta:x_end,y_sta:y_end,z_sta:z_end,1:ns))
  ! allocate(sf_avg(x_sta:x_end,y_sta:y_end,z_sta:z_end,1:3))
  ! allocate(ri(x_sta:x_end,y_sta:y_end,z_sta:z_end))
  ! allocate(rayleighindex(x_sta:x_end,y_sta:y_end,z_sta:z_end))
  !allocate(q(1-ibd:nx+ibd,1-jbd:ny+jbd,1-kbd:nz+kbd,4))

  return
end subroutine

subroutine deallocation
  use globals

  deallocate(u_local)
  deallocate(v_local)
  deallocate(w_local)
  deallocate(r_local)
  deallocate(p_local)
  deallocate(t_local)
  deallocate(h_local)
  deallocate(y_local)
 
end subroutine deallocation

subroutine Z_ini
  use globals
  implicit none
  logical::oxi_flag
  double precision :: xn2
  
  WC  = 12.0112d0
  WO  = 15.9994d0
  WH  = 1.00797d0
  WN  = 14.0067d0

  owc=1d0/wc
  owh=1d0/wh
  own=1d0/wn
  owo=1d0/wo

  open(111,file="./dat/z_setting.dat",form='formatted')
  read(111,*) cfuel
  read(111,*) hfuel
  read(111,*) nfuel
  read(111,*) ofuel
  read(111,*) zc1, zc2
  read(111,*) zh1, zh2
  read(111,*) zn1, zn2
  read(111,*) zo1, zo2
  read(111,*) xo2
  read(111,*) oxi_flag
  read(111,*) zst
  close(111)
  write(*,*) "ZC1, Zc2",zc1, zc2
  write(*,*) "Zh1, Zh2",zh1, zh2
  write(*,*) "Zn1, Zn2",zn1, zn2
  write(*,*) "Zo1, Zo2",zo1, zo2

  if(oxi_flag)then
     xn2=1d0-xo2
     zo2=xo2*WO / (xo2*WO+ xn2*WN)
     zn2=xn2*WN/ (xo2*WO+ xn2*WN)
  endif

  nu_o2= cfuel + hfuel*0.25d0 - ofuel*0.5d0

  if(cfuel/=0d0)then
     owc=owc/cfuel
  else
     owc=0d0
  endif
  if(hfuel/=0d0)then
     owh=owh/hfuel
  else
     owh=0d0
  endif
  if(nfuel/=0d0)then
     own=own/nfuel
  else
     own=0d0
  endif
  owo=owo/nu_o2
  
  allocate(coeffZ(nf,4))

  open(222,file='./dat/coeff4calcZ.dat',form='formatted')
  read(222,'()')
  do i=1,nf
   read(222,*)  coeffZ(i,1:4) ! C, H, N, O
   if(i==nf-5) write(*,*) "coeffZ",coeffZ(i,1:4)
  enddo
  close(222)

  write(*,*) "coeff_Z_H"
  write(*,*) coeffZ(:,2)
   
  
  return
end subroutine Z_ini

subroutine calc_z(ysub,zsub,zhsub,znsub,hn_ratio)
  use globals 
  double precision,intent(in)::ysub(1:nf)
  double precision,intent(out)::zsub,zhsub,znsub,hn_ratio

  double precision:: Zn_fuel
  
  zc=0d0
  zh=0d0
  zo=0d0
  zn=0d0
  
  do ii=1,nf
     Zc=Zc+ysub(ii)*coeffZ(ii,1)
     Zh=Zh+ysub(ii)*coeffZ(ii,2)
     Zn=Zn+ysub(ii)*coeffZ(ii,3)
     Zo=Zo+ysub(ii)*coeffZ(ii,4)
  enddo

  zsub = ((Zc -Zc2)*owc+(Zh -Zh2)*owh+(Zo2-Zo )*owo + (Zn-Zn2)*own) / &
      ((Zc1-Zc2)*owc+(Zh1-Zh2)*owh+(Zo2-Zo1)*owo + (Zn1-Zn2)*own)
  
  zhsub = (Zh-Zh2)/(Zh1-Zh2)
  if(zhsub<0d0)then
     write(*,*) "Zh<0d0"
     write(*,*) "  Zhsub",zhsub
     write(*,*) "     Zh",zh
     write(*,*) "    Zh1",zh1
     write(*,*) "    Zh2",zh2
     write(*,*) "    ysub",ysub
  endif
  
  znsub = (Zn-Zn2)/(Zn1-Zn2)
  
  Zn_fuel = Zn - Zo * (1d0-xo2)*wn*2d0/(xo2*wo*2d0)

  if(zn_fuel>1d-14)then
     hn_ratio= (zh/wh) / (zn_fuel/wn)
  else
     hn_ratio = 0d0
  endif
  return
end subroutine calc_z

subroutine chem_ini
  use globals
  implicit none
  integer:: njj
  ifn=111

  allocate(wg(1:nf))
  allocate(cmu(1:nf,1:4))
  allocate (clam(1:nf,1:4))
  allocate (cdif(1:nf,1:nf,1:4))
  allocate (ccpl(1:nf,1:7))
  allocate (ccph(1:nf,1:7))
  allocate (tmpr(1:nf,1:3))
  
  
  open(ifn,file=trim(adjustl(chem_dir))//'ck.out',form='formatted')
  read(ifn,*) njj
  read(ifn,*) wg
  read(ifn,*) cmu
  read(ifn,*) clam
  read(ifn,*) cdif
  read(ifn,*) species_names_fk3(1:njj)
  close(ifn)

  open(ifn,file=trim(adjustl(chem_dir))//'therm.out',form='formatted')
  read(ifn,*) njj
  read(ifn,*) tmpr
  read(ifn,*) ccpl,ccph
  close(ifn)

  ugconst=8.314472d0
    
  return
end subroutine chem_ini

subroutine calc_Dh(y,r,temp,dh)
  use globals, only : nf,wg,clam,tmpr,ccpl,ccph,ugconst
  implicit none
  double precision, intent(in)::y(1:nf),r,temp
  double precision, intent(out)::dh
  double precision:: cpa,cp, lam0(1:nf),lam, lnt,x_chem(1:nf),dum1,dum2
  integer:: i,l
  
  cp=0d0
  do l=1,nf
    if(temp>=tmpr(l,3))then
       cpa = ccph(l,5)*temp**4 + ccph(l,4)*temp**3 + ccph(l,3)*temp**2 &
            +ccph(l,2)*temp**1 + ccph(l,1)
    else
       cpa = ccpl(l,5)*temp**4 + ccpl(l,4)*temp**3 + ccpl(l,3)*temp**2 &
            +ccpl(l,2)*temp**1 + ccpl(l,1)       
    end if
    cp=cp+cpa*ugconst/wg(l)*y(l)*1000.d0
  enddo

  dum1=sum(y/wg)
  do i=1,nf
     x_chem(i)=y(i)/(wg(i)*dum1)
  enddo

  lam=0.d0
  lam0=0.d0
  lnt=dlog(temp)
  do i=1,nf
     lam0(i)=clam(i,1)+lnt*(clam(i,2)+lnt*(clam(i,3)+lnt*clam(i,4)))
  enddo
  lam0=dexp(lam0)

 dum1=sum(x_chem*lam0)
 dum2=sum(x_chem/lam0)
 lam=(dum1+1.d0/dum2)*0.5d0/1.d5 
  

  dh = lam/cp/r
  
  return
end subroutine calc_Dh
