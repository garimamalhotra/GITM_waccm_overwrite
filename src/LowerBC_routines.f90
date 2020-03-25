subroutine load_files_mem_eff
    use netcdf
    use ModTime, ONLY: CurrentTime
    use ModInputs
    use ModLowerBC
    
    implicit none
    
    integer :: alt_dimid,lat_dimid,lon_dimid,hour_dimid,ncid,hour_id
    !To enquire about the netCDF to know more about an unknown netCDF
    integer ::ndims_in,nvars_in,ngatts_in,unlimdimid_in,ilow,ihigh,real_hours
    integer :: iError,nhours
    integer :: i,l=0
    
    real_hours=24*nBCFiles
   
    if (.not. allocated(hour_arr)) then
      do i=1,nBCFiles       
          !Open the netcdf file to read the data  
          iError=nf90_open(BCFileNames(i),nf90_nowrite,ncid)
          iError=nf90_inquire(ncid,ndims_in,nvars_in,ngatts_in,unlimdimid_in) 
          ! Getting the dimensions
          iError=nf90_inq_dimid(ncid,"nhours",hour_dimid)
          iError=nf90_Inquire_Dimension(ncid,hour_dimid,len=nHours)

          ilow=((i-1)*nHours)+1
          ihigh=(i*nHours)

          if (i==1) then
            allocate(hour_arr(real_hours))
          endif

          iError=nf90_inq_varid(ncid,"Time",hour_id)
          iError=nf90_get_var(ncid,hour_id,hour_arr(ilow:ihigh))
      enddo
    endif

    do l=1,real_hours
        if (hour_arr(l) > CurrentTime) then
            exit
        end if
    end do
    !deallocate(hour_arr) 
    !write(*,*) 'Hour Arr : ', hour_arr
    !write(*,*) 'Final l ', l, hour_arr(l-1), CurrentTime, hour_arr(l) 
    !write(*,*) 'Stepping hour : ', hour_arr(1),&
    !hour_arr(24),hour_arr(25),hour_arr(48),&
    !hour_arr(49),hour_arr(72), hour_arr(73)
    call load_file_w_indices_only_one(l)

end subroutine load_files_mem_eff

subroutine load_file_w_indices_only_one(l)
    use ModLowerBC
    use ModInputs
    use netcdf
    use ModTime, ONLY: CurrentTime

    implicit none
    integer, intent(in) ::l
    integer :: file_index, index_within_file
    integer :: i=0
    integer :: nLats,nLons,nHours,num_Alts  ! Changed nAlts to num_Alts because of
    !namespace conflict with ModInput -> ModPlanet.f90
    real,dimension(24) :: temp
    integer :: ncid,alt_id,lat_id,lon_id,hour_id,&
               ktemp_id,zonal_id,merid_id,density_id,o_id,o2_id,n2_id,&
               co2_id,no_id,n_id,vert_id,iError
    integer :: alt_dimid,lat_dimid,lon_dimid,hour_dimid
    !To enquire about the netCDF to know more about an unknown netCDF
    integer::ndims_in,nvars_in,ngatts_in,unlimdimid_in,e_Error

    !To enquire about the netCDF to know more about an unknown netCDF
    integer :: real_hours,k


    !This 'l' is the index of the next element - between which the current time lies
    file_index=(l/24)+1 !because we don't want 0 index for the first file
    index_within_file=modulo(l,24)
    !write(*,*) 'value of l index : ' ,l, file_index, index_within_file
    k=0
    real_hours=2

    if (index_within_file == 1) then
        do i=file_index-1,file_index
          if (.not. allocated(hour)) then
              call allocate_vars(i,ncid)
          end if

          !write(*,*) 'before time : ',CurrentTime,allocated(hour), hour(1), hour(real_hours)
          if (CurrentTime >= hour(real_hours)) then
              iError=nf90_open(BCFileNames(i),nf90_nowrite,ncid)
              iError=nf90_inquire(ncid,ndims_in,nvars_in,ngatts_in,unlimdimid_in) 
              ! Getting the dimensions
              iError=nf90_inq_dimid(ncid,"nlats",lat_dimid)
              iError=nf90_inq_dimid(ncid,"nlons",lon_dimid)
              iError=nf90_inq_dimid(ncid,"nalts",alt_dimid)
              iError=nf90_inq_dimid(ncid,"nhours",hour_dimid)
              iError=nf90_Inquire_Dimension(ncid,lat_dimid,len=nLats)
              iError=nf90_Inquire_Dimension(ncid,lon_dimid,len=nLons)
              iError=nf90_Inquire_Dimension(ncid,alt_dimid,len=num_Alts)
              iError=nf90_Inquire_Dimension(ncid,hour_dimid,len=nHours)

              !write(*,*) 'allocating : ', allocated(hour), 'time : ',CurrentTime, hour ,i,l
              
              iError=nf90_inq_varid(ncid,"Latitude",lat_id)
              iError=nf90_inq_varid(ncid,"Longitude",lon_id)
              iError=nf90_inq_varid(ncid,"Altitude",alt_id)
              iError=nf90_inq_varid(ncid,"Time",hour_id)
              iError=nf90_inq_varid(ncid,"Temperature",ktemp_id)
              iError=nf90_inq_varid(ncid,"Zonal",zonal_id)
              iError=nf90_inq_varid(ncid,"Meridional",merid_id)
              iError=nf90_inq_varid(ncid,"Density",density_id)
              iError=nf90_inq_varid(ncid,"O",o_id)
              iError=nf90_inq_varid(ncid,"O2",o2_id)
              iError=nf90_inq_varid(ncid,"N2",n2_id)
              iError=nf90_inq_varid(ncid,"CO2",co2_id)
              iError=nf90_inq_varid(ncid,"NO",no_id)
              iError=nf90_inq_varid(ncid,"N",n_id)
              iError=nf90_inq_varid(ncid,"Omega",vert_id)
              

              if (k==0) then
                  !iError=nf90_get_var(ncid,lat_id,lat_waccm)
                  !iError=nf90_get_var(ncid,lon_id,lon_waccm)
                  !iError=nf90_get_var(ncid,alt_id,alt_waccm)
                  iError=nf90_get_var(ncid,hour_id,temp)
                  hour(1)=temp(nHours)
                  iError=nf90_get_var(ncid,ktemp_id,ktemp(:,:,:,1),start=(/ 1,1,1,nHours /),count=(/ nLons,nLats,num_Alts,1 /))
                  iError=nf90_get_var(ncid,zonal_id,zonal(:,:,:,1),start=(/ 1,1,1,nHours /),count=(/ nLons,nLats,num_Alts,1 /))
                  iError=nf90_get_var(ncid,merid_id,merid(:,:,:,1),start=(/ 1,1,1,nHours /),count=(/ nLons,nLats,num_Alts,1 /))
                  iError=nf90_get_var(ncid,density_id,neutral_density(:,:,:,1),start=(/ 1,1,1,nHours /),count=(/ nLons,nLats,num_Alts,1 /))
                  iError=nf90_get_var(ncid,o_id,o(:,:,:,1),start=(/ 1,1,1,nHours /),count=(/ nLons,nLats,num_Alts,1 /))
                  iError=nf90_get_var(ncid,o2_id,o2(:,:,:,1),start=(/ 1,1,1,nHours /),count=(/ nLons,nLats,num_Alts,1 /))
                  iError=nf90_get_var(ncid,n2_id,nitrogen(:,:,:,1),start=(/ 1,1,1,nHours /),count=(/ nLons,nLats,num_Alts,1 /))
                  iError=nf90_get_var(ncid,co2_id,co2(:,:,:,1),start=(/ 1,1,1,nHours /),count=(/ nLons,nLats,num_Alts,1 /))
                  iError=nf90_get_var(ncid,no_id,no(:,:,:,1),start=(/ 1,1,1,nHours /),count=(/ nLons,nLats,num_Alts,1 /))
                  iError=nf90_get_var(ncid,n_id,atomic_n(:,:,:,1),start=(/ 1,1,1,nHours /),count=(/ nLons,nLats,num_Alts,1 /))
                  iError=nf90_get_var(ncid,vert_id,vert(:,:,:,1),start=(/ 1,1,1,nHours /),count=(/ nLons,nLats,num_Alts,1 /))
                  !write(*,*) ,'k = ', k, shape(ktemp), hour(1), CurrentTime,'index : ', l,file_index, index_within_file 
              else
                  iError=nf90_get_var(ncid,lat_id,lat_waccm)
                  iError=nf90_get_var(ncid,lon_id,lon_waccm)
                  iError=nf90_get_var(ncid,alt_id,alt_waccm)
                  iError=nf90_get_var(ncid,hour_id,temp)
                  hour(2)=temp(1)
                  iError=nf90_get_var(ncid,ktemp_id,ktemp(:,:,:,2),start=(/ 1,1,1,1 /),count=(/ nLons,nLats,num_Alts,1 /))
                  iError=nf90_get_var(ncid,zonal_id,zonal(:,:,:,2),start=(/ 1,1,1,1 /),count=(/ nLons,nLats,num_Alts,1 /))
                  iError=nf90_get_var(ncid,merid_id,merid(:,:,:,2),start=(/ 1,1,1,1 /),count=(/ nLons,nLats,num_Alts,1 /))
                  iError=nf90_get_var(ncid,density_id,neutral_density(:,:,:,2),start=(/ 1,1,1,1 /),count=(/ nLons,nLats,num_Alts,1 /))
                  iError=nf90_get_var(ncid,o_id,o(:,:,:,2),start=(/ 1,1,1,1 /),count=(/ nLons,nLats,num_Alts,1 /))
                  iError=nf90_get_var(ncid,o2_id,o2(:,:,:,2),start=(/ 1,1,1,1 /),count=(/ nLons,nLats,num_Alts,1 /))
                  iError=nf90_get_var(ncid,n2_id,nitrogen(:,:,:,2),start=(/ 1,1,1,1 /),count=(/ nLons,nLats,num_Alts,1 /))
                  iError=nf90_get_var(ncid,co2_id,co2(:,:,:,2),start=(/ 1,1,1,1 /),count=(/ nLons,nLats,num_Alts,1 /))
                  iError=nf90_get_var(ncid,no_id,no(:,:,:,2),start=(/ 1,1,1,1 /),count=(/ nLons,nLats,num_Alts,1 /))
                  iError=nf90_get_var(ncid,n_id,atomic_n(:,:,:,2),start=(/ 1,1,1,1 /),count=(/ nLons,nLats,num_Alts,1 /))
                  iError=nf90_get_var(ncid,vert_id,vert(:,:,:,2),start=(/ 1,1,1,1 /),count=(/ nLons,nLats,num_Alts,1 /))
                  !write(*,*) ,'k = ', k, shape(ktemp), hour(2),CurrentTime,'index : ', l,file_index, index_within_file 


              end if
              !write(*,*) 'after time : ', k, 'hour, ',hour,shape(hour),i,l,'index : ', file_index, index_within_file
              !write(*,*) 'starting of file : ', hour(1), CurrentTime, hour(2)
              iError=nf90_close(ncid)
          end if

          k=k+1 
        end do
    else
        i=file_index
        if (.not. allocated(hour)) then
          call allocate_vars(i,ncid)
        end if

        if (CurrentTime >= hour(real_hours)) then
          iError=nf90_open(BCFileNames(i),nf90_nowrite,ncid)
          iError=nf90_inquire(ncid,ndims_in,nvars_in,ngatts_in,unlimdimid_in) 
          ! Getting the dimensions
          iError=nf90_inq_dimid(ncid,"nlats",lat_dimid)
          iError=nf90_inq_dimid(ncid,"nlons",lon_dimid)
          iError=nf90_inq_dimid(ncid,"nalts",alt_dimid)
          iError=nf90_inq_dimid(ncid,"nhours",hour_dimid)
          iError=nf90_Inquire_Dimension(ncid,lat_dimid,len=nLats)
          iError=nf90_Inquire_Dimension(ncid,lon_dimid,len=nLons)
          iError=nf90_Inquire_Dimension(ncid,alt_dimid,len=num_Alts)
          iError=nf90_Inquire_Dimension(ncid,hour_dimid,len=nHours)

          iError=nf90_inq_varid(ncid,"Latitude",lat_id)
          iError=nf90_inq_varid(ncid,"Longitude",lon_id)
          iError=nf90_inq_varid(ncid,"Altitude",alt_id)
          iError=nf90_inq_varid(ncid,"Time",hour_id)
          iError=nf90_inq_varid(ncid,"Temperature",ktemp_id)
          iError=nf90_inq_varid(ncid,"Zonal",zonal_id)
          iError=nf90_inq_varid(ncid,"Meridional",merid_id)
          iError=nf90_inq_varid(ncid,"Density",density_id)
          iError=nf90_inq_varid(ncid,"O",o_id)
          iError=nf90_inq_varid(ncid,"O2",o2_id)
          iError=nf90_inq_varid(ncid,"N2",n2_id)
          iError=nf90_inq_varid(ncid,"CO2",co2_id)
          iError=nf90_inq_varid(ncid,"NO",no_id)
          iError=nf90_inq_varid(ncid,"N",n_id)
          iError=nf90_inq_varid(ncid,"Omega",vert_id)
 
          iError=nf90_get_var(ncid,lat_id,lat_waccm)
          iError=nf90_get_var(ncid,lon_id,lon_waccm)
          iError=nf90_get_var(ncid,alt_id,alt_waccm)
          iError=nf90_get_var(ncid,hour_id,temp)
          hour(1:2)=temp(index_within_file-1:index_within_file)
          !write(*,*),'Before assign : ', ktemp
          iError=nf90_get_var(ncid,ktemp_id,ktemp,start=(/ 1,1,1,index_within_file-1 /),count=(/ nLons,nLats,num_Alts,2 /))
          iError=nf90_get_var(ncid,zonal_id,zonal,start=(/ 1,1,1,index_within_file-1 /),count=(/ nLons,nLats,num_Alts,2 /))
          iError=nf90_get_var(ncid,merid_id,merid,start=(/ 1,1,1,index_within_file-1 /),count=(/ nLons,nLats,num_Alts,2 /))
          iError=nf90_get_var(ncid,density_id,neutral_density,start=(/ 1,1,1,index_within_file-1 /),count=(/ nLons,nLats,num_Alts,2 /))
          iError=nf90_get_var(ncid,o_id,o,start=(/ 1,1,1,index_within_file-1 /),count=(/ nLons,nLats,num_Alts,2 /))
          iError=nf90_get_var(ncid,o2_id,o2,start=(/ 1,1,1,index_within_file-1 /),count=(/ nLons,nLats,num_Alts,2 /))
          iError=nf90_get_var(ncid,n2_id,nitrogen,start=(/ 1,1,1,index_within_file-1 /),count=(/ nLons,nLats,num_Alts,2 /))
          iError=nf90_get_var(ncid,co2_id,co2,start=(/ 1,1,1,index_within_file-1 /),count=(/ nLons,nLats,num_Alts,2 /))
          iError=nf90_get_var(ncid,no_id,no,start=(/ 1,1,1,index_within_file-1 /),count=(/ nLons,nLats,num_Alts,2 /))
          iError=nf90_get_var(ncid,n_id,atomic_n,start=(/ 1,1,1,index_within_file-1 /),count=(/ nLons,nLats,num_Alts,2 /))
          iError=nf90_get_var(ncid,vert_id,vert,start=(/ 1,1,1,index_within_file-1 /),count=(/ nLons,nLats,num_Alts,2 /))
          !write(*,*) ,'In between', hour(1), CurrentTime,hour(2),'index : ', index_within_file, file_index, l
          !write(*,*) 'After assign : ', ktemp
        end if
    end if
    
            
end subroutine load_file_w_indices_only_one

subroutine allocate_vars(i)
      use ModLowerBC
      use ModInputs
      use netcdf
      ! Name of the dimensionsnum_Alts,
      implicit none

      integer :: ncid,iError,alt_dimid,lat_dimid,lon_dimid,hour_dimid
      !To enquire about the netCDF to know more about an unknown netCDF
      integer::ndims_in,nvars_in,ngatts_in,unlimdimid_in,e_Error,real_hours
      integer :: nLats,nLons,nHours,num_Alts 
      integer, intent(in) ::  i

      real_hours=2 

      iError=nf90_open(BCFileNames(i),nf90_nowrite,ncid)
      iError=nf90_inquire(ncid,ndims_in,nvars_in,ngatts_in,unlimdimid_in) 
      ! Getting the dimensions
      iError=nf90_inq_dimid(ncid,"nlats",lat_dimid)
      iError=nf90_inq_dimid(ncid,"nlons",lon_dimid)
      iError=nf90_inq_dimid(ncid,"nalts",alt_dimid)
      iError=nf90_inq_dimid(ncid,"nhours",hour_dimid)
      iError=nf90_Inquire_Dimension(ncid,lat_dimid,len=nLats)
      iError=nf90_Inquire_Dimension(ncid,lon_dimid,len=nLons)
      iError=nf90_Inquire_Dimension(ncid,alt_dimid,len=num_Alts)
      iError=nf90_Inquire_Dimension(ncid,hour_dimid,len=nHours)
    
      !ilow=((k-1)*nHours)+1
      !ihigh=(k*nHours)
      !write(*,*) 'step 3', ilow, ihigh
      !real_hours=(nHours*nBCFiles_indices)+1
      ! For variables and data
      if (.not. allocated(hour)) then
          allocate(lat_waccm(nLats))
          allocate(lon_waccm(nLons))
          allocate(alt_waccm(num_Alts))
          allocate(hour(real_hours))
          allocate(ktemp(nLons,nLats,num_Alts,real_hours))
          allocate(zonal(nLons,nLats,num_Alts,real_hours))
          allocate(merid(nLons,nLats,num_Alts,real_hours))
          allocate(neutral_density(nLons,nLats,num_Alts,real_hours))
          allocate(o(nLons,nLats,num_Alts,real_hours))
          allocate(o2(nLons,nLats,num_Alts,real_hours))
          allocate(nitrogen(nLons,nLats,num_Alts,real_hours))
          allocate(co2(nLons,nLats,num_Alts,real_hours))
          allocate(no(nLons,nLats,num_Alts,real_hours))
          allocate(atomic_n(nLons,nLats,num_Alts,real_hours))
          allocate(vert(nLons,nLats,num_Alts,real_hours))
        
          allocate(interp_ktemp(nLons,nLats,num_Alts))
          allocate(interp_zonal(nLons,nLats,num_Alts))
          allocate(interp_merid(nLons,nLats,num_Alts))
          allocate(interp_density(nLons,nLats,num_Alts))
          allocate(interp_o(nLons,nLats,num_Alts))
          allocate(interp_o2(nLons,nLats,num_Alts))
          allocate(interp_nitrogen(nLons,nLats,num_Alts))
          allocate(interp_co2(nLons,nLats,num_Alts))
          allocate(interp_no(nLons,nLats,num_Alts))
          allocate(interp_atomic_n(nLons,nLats,num_Alts))
          allocate(interp_vert(nLons,nLats,num_Alts))
          hour=0.0
      end if
      !write(*,*) 'Allocated : ', BCFileNames(i)

end subroutine allocate_vars

subroutine load_file_w_indices(l)
    use ModLowerBC
    use ModInputs
    use netcdf
    use ModTime, ONLY: CurrentTime

    implicit none
    integer, intent(in) ::l
    integer :: file_index
    integer :: i=0
    integer :: nLats,nLons,nHours,num_Alts  ! Changed nAlts to num_Alts because of
    !namespace conflict with ModInput -> ModPlanet.f90

    ! Name of the dimensionsnum_Alts,
    integer :: ncid,iError,alt_id,lat_id,lon_id,hour_id
    integer :: ktemp_id,zonal_id,merid_id,density_id,o_id,o2_id,n2_id,&
          co2_id,no_id,n_id,vert_id
    integer :: alt_dimid,lat_dimid,lon_dimid,hour_dimid
    !To enquire about the netCDF to know more about an unknown netCDF
    integer::ndims_in,nvars_in,ngatts_in,unlimdimid_in,e_Error,ilow,ihigh,real_hours,k
    integer :: nBCFiles_indices=1

    file_index=(l/24)+1
    !write(*,*) 'value of l index : ' ,l, file_index, BCFilenames(file_index)
    !write(*,*) 'step 1 ', file_index
    k=1
    do i=file_index-1,file_index 
          iError=nf90_open(BCFileNames(i),nf90_nowrite,ncid)
          iError=nf90_inquire(ncid,ndims_in,nvars_in,ngatts_in,unlimdimid_in) 
          ! Getting the dimensions
          iError=nf90_inq_dimid(ncid,"nlats",lat_dimid)
          iError=nf90_inq_dimid(ncid,"nlons",lon_dimid)
          iError=nf90_inq_dimid(ncid,"nalts",alt_dimid)
          iError=nf90_inq_dimid(ncid,"nhours",hour_dimid)
          iError=nf90_Inquire_Dimension(ncid,lat_dimid,len=nLats)
          iError=nf90_Inquire_Dimension(ncid,lon_dimid,len=nLons)
          iError=nf90_Inquire_Dimension(ncid,alt_dimid,len=num_Alts)
          iError=nf90_Inquire_Dimension(ncid,hour_dimid,len=nHours)
        
          ilow=((k-1)*nHours)+1
          ihigh=(k*nHours)
          !write(*,*) 'step 3', ilow, ihigh
          real_hours=nHours*nBCFiles_indices

          ! For variables and data
          if (.not. allocated(hour)) then
              allocate(lat_waccm(nLats))
              allocate(lon_waccm(nLons))
              allocate(alt_waccm(num_Alts))
              allocate(hour(real_hours))
              allocate(ktemp(nLons,nLats,num_Alts,real_hours))
              allocate(zonal(nLons,nLats,num_Alts,real_hours))
              allocate(merid(nLons,nLats,num_Alts,real_hours))
              allocate(neutral_density(nLons,nLats,num_Alts,real_hours))
              allocate(o(nLons,nLats,num_Alts,real_hours))
              allocate(o2(nLons,nLats,num_Alts,real_hours))
              allocate(nitrogen(nLons,nLats,num_Alts,real_hours))
              allocate(co2(nLons,nLats,num_Alts,real_hours))
              allocate(no(nLons,nLats,num_Alts,real_hours))
              allocate(atomic_n(nLons,nLats,num_Alts,real_hours))
              allocate(vert(nLons,nLats,num_Alts,real_hours))
            
              allocate(interp_ktemp(nLons,nLats,num_Alts))
              allocate(interp_zonal(nLons,nLats,num_Alts))
              allocate(interp_merid(nLons,nLats,num_Alts))
              allocate(interp_density(nLons,nLats,num_Alts))
              allocate(interp_o(nLons,nLats,num_Alts))
              allocate(interp_o2(nLons,nLats,num_Alts))
              allocate(interp_nitrogen(nLons,nLats,num_Alts))
              allocate(interp_co2(nLons,nLats,num_Alts))
              allocate(interp_no(nLons,nLats,num_Alts))
              allocate(interp_atomic_n(nLons,nLats,num_Alts))
              allocate(interp_vert(nLons,nLats,num_Alts))
              hour=0.0
          end if
          !write(*,*) 'before time : ',CurrentTime,allocated(hour), hour(1), hour(real_hours)
          if (CurrentTime >= hour(real_hours)) then 
              write(*,*) 'allocated : ', allocated(hour), 'time : ',CurrentTime, hour(real_hours),i,l
              iError=nf90_inq_varid(ncid,"Latitude",lat_id)
              iError=nf90_inq_varid(ncid,"Longitude",lon_id)
              iError=nf90_inq_varid(ncid,"Altitude",alt_id)
              iError=nf90_inq_varid(ncid,"Time",hour_id)
              iError=nf90_inq_varid(ncid,"Temperature",ktemp_id)
              iError=nf90_inq_varid(ncid,"Zonal",zonal_id)
              iError=nf90_inq_varid(ncid,"Meridional",merid_id)
              iError=nf90_inq_varid(ncid,"Density",density_id)
              iError=nf90_inq_varid(ncid,"O",o_id)
              iError=nf90_inq_varid(ncid,"O2",o2_id)
              iError=nf90_inq_varid(ncid,"N2",n2_id)
              iError=nf90_inq_varid(ncid,"CO2",co2_id)
              iError=nf90_inq_varid(ncid,"NO",no_id)
              iError=nf90_inq_varid(ncid,"N",n_id)
              iError=nf90_inq_varid(ncid,"Omega",vert_id)

              iError=nf90_get_var(ncid,lat_id,lat_waccm)
              iError=nf90_get_var(ncid,lon_id,lon_waccm)
              iError=nf90_get_var(ncid,alt_id,alt_waccm)
              iError=nf90_get_var(ncid,hour_id,hour(ilow:ihigh))
              iError=nf90_get_var(ncid,ktemp_id,ktemp(:,:,:,ilow:ihigh))
              iError=nf90_get_var(ncid,zonal_id,zonal(:,:,:,ilow:ihigh))
              iError=nf90_get_var(ncid,merid_id,merid(:,:,:,ilow:ihigh))
              iError=nf90_get_var(ncid,density_id,neutral_density(:,:,:,ilow:ihigh))
              iError=nf90_get_var(ncid,o_id,o(:,:,:,ilow:ihigh))
              iError=nf90_get_var(ncid,o2_id,o2(:,:,:,ilow:ihigh))
              iError=nf90_get_var(ncid,n2_id,nitrogen(:,:,:,ilow:ihigh))
              iError=nf90_get_var(ncid,co2_id,co2(:,:,:,ilow:ihigh))
              iError=nf90_get_var(ncid,no_id,no(:,:,:,ilow:ihigh))
              iError=nf90_get_var(ncid,n_id,atomic_n(:,:,:,ilow:ihigh))
              iError=nf90_get_var(ncid,vert_id,vert(:,:,:,ilow:ihigh))
              write(*,*) 'after time ilow : ', ilow, 'ihigh',ihigh, real_hours,shape(hour),'final hour', hour(real_hours) 
              iError=nf90_close(ncid)
          endif

          k=k+1 
    end do        
end subroutine load_file_w_indices

subroutine read_netcdf
  use netcdf
  use ModLowerBC
  use ModInputs
  
  implicit none
  integer :: i=0
  integer :: nLats,nLons,nHours,num_Alts  ! Changed nAlts to num_Alts because of
  !namespace conflict with ModInput -> ModPlanet.f90

  ! Name of the dimensionsnum_Alts,
  integer :: ncid,iError,alt_id,lat_id,lon_id,hour_id
  integer :: ktemp_id,zonal_id,merid_id,density_id,o_id,o2_id,n2_id,&
      co2_id,no_id,n_id,vert_id
  integer :: alt_dimid,lat_dimid,lon_dimid,hour_dimid
  !To enquire about the netCDF to know more about an unknown netCDF
  integer::ndims_in,nvars_in,ngatts_in,unlimdimid_in,e_Error,ilow,ihigh,real_hours
  integer,dimension(1):: shape_arr

  do i=1,nBCFiles         
      !Open the netcdf file to read the data  
      iError=nf90_open(BCFileNames(i),nf90_nowrite,ncid)
      iError=nf90_inquire(ncid,ndims_in,nvars_in,ngatts_in,unlimdimid_in) 
      ! Getting the dimensions
      iError=nf90_inq_dimid(ncid,"nlats",lat_dimid)
      iError=nf90_inq_dimid(ncid,"nlons",lon_dimid)
      iError=nf90_inq_dimid(ncid,"nalts",alt_dimid)
      iError=nf90_inq_dimid(ncid,"nhours",hour_dimid)
      iError=nf90_Inquire_Dimension(ncid,lat_dimid,len=nLats)
      iError=nf90_Inquire_Dimension(ncid,lon_dimid,len=nLons)
      iError=nf90_Inquire_Dimension(ncid,alt_dimid,len=num_Alts)
      iError=nf90_Inquire_Dimension(ncid,hour_dimid,len=nHours)
    
      ilow=((i-1)*nHours)+1
      ihigh=(i*nHours)

      ! For variables and data
      if (i==1) then
          real_hours=nHours*nBCFiles

          allocate(lat_waccm(nLats))
          allocate(lon_waccm(nLons))
          allocate(alt_waccm(num_Alts))
          allocate(hour(real_hours))
          allocate(ktemp(nLons,nLats,num_Alts,real_hours))
          allocate(zonal(nLons,nLats,num_Alts,real_hours))
          allocate(merid(nLons,nLats,num_Alts,real_hours))
          allocate(neutral_density(nLons,nLats,num_Alts,real_hours))
          allocate(o(nLons,nLats,num_Alts,real_hours))
          allocate(o2(nLons,nLats,num_Alts,real_hours))
          allocate(nitrogen(nLons,nLats,num_Alts,real_hours))
          allocate(co2(nLons,nLats,num_Alts,real_hours))
          allocate(no(nLons,nLats,num_Alts,real_hours))
          allocate(atomic_n(nLons,nLats,num_Alts,real_hours))
          allocate(vert(nLons,nLats,num_Alts,real_hours))
        
          
          allocate(interp_ktemp(nLons,nLats,num_Alts))
          allocate(interp_zonal(nLons,nLats,num_Alts))
          allocate(interp_merid(nLons,nLats,num_Alts))
          allocate(interp_density(nLons,nLats,num_Alts))
          allocate(interp_o(nLons,nLats,num_Alts))
          allocate(interp_o2(nLons,nLats,num_Alts))
          allocate(interp_nitrogen(nLons,nLats,num_Alts))
          allocate(interp_co2(nLons,nLats,num_Alts))
          allocate(interp_no(nLons,nLats,num_Alts))
          allocate(interp_atomic_n(nLons,nLats,num_Alts))
          allocate(interp_vert(nLons,nLats,num_Alts))

      end if

      iError=nf90_inq_varid(ncid,"Latitude",lat_id)
      iError=nf90_inq_varid(ncid,"Longitude",lon_id)
      iError=nf90_inq_varid(ncid,"Altitude",alt_id)
      iError=nf90_inq_varid(ncid,"Time",hour_id)
      iError=nf90_inq_varid(ncid,"Temperature",ktemp_id)
      iError=nf90_inq_varid(ncid,"Zonal",zonal_id)
      iError=nf90_inq_varid(ncid,"Meridional",merid_id)
      iError=nf90_inq_varid(ncid,"Density",density_id)
      iError=nf90_inq_varid(ncid,"O",o_id)
      iError=nf90_inq_varid(ncid,"O2",o2_id)
      iError=nf90_inq_varid(ncid,"N2",n2_id)
      iError=nf90_inq_varid(ncid,"CO2",co2_id)
      iError=nf90_inq_varid(ncid,"NO",no_id)
      iError=nf90_inq_varid(ncid,"N",n_id)
      iError=nf90_inq_varid(ncid,"Vert",vert_id)

      iError=nf90_get_var(ncid,lat_id,lat_waccm)
      iError=nf90_get_var(ncid,lon_id,lon_waccm)
      iError=nf90_get_var(ncid,alt_id,alt_waccm)
      iError=nf90_get_var(ncid,hour_id,hour(ilow:ihigh))
      iError=nf90_get_var(ncid,ktemp_id,ktemp(:,:,:,ilow:ihigh))
      iError=nf90_get_var(ncid,zonal_id,zonal(:,:,:,ilow:ihigh))
      iError=nf90_get_var(ncid,merid_id,merid(:,:,:,ilow:ihigh))
      iError=nf90_get_var(ncid,density_id,neutral_density(:,:,:,ilow:ihigh))
      iError=nf90_get_var(ncid,o_id,o(:,:,:,ilow:ihigh))
      iError=nf90_get_var(ncid,o2_id,o2(:,:,:,ilow:ihigh))
      iError=nf90_get_var(ncid,n2_id,nitrogen(:,:,:,ilow:ihigh))
      iError=nf90_get_var(ncid,co2_id,co2(:,:,:,ilow:ihigh))
      iError=nf90_get_var(ncid,no_id,no(:,:,:,ilow:ihigh))
      iError=nf90_get_var(ncid,n_id,atomic_n(:,:,:,ilow:ihigh))
      iError=nf90_get_var(ncid,vert_id,vert(:,:,:,ilow:ihigh))
   
      iError=nf90_close(ncid)

     
  end do         
end subroutine read_netcdf 

subroutine interp_3Dspace(Alt,Lat,Lon,interp_ktemp_space,interp_zonal_space,&
        interp_merid_space,interp_density_space,interp_o_space,interp_o2_space,&
        interp_nitrogen_space,interp_co2_space,interp_no_space,interp_n_space,interp_vert_space)
        
    use ModLowerBC
    
    implicit none

    real, intent(in) :: Alt,Lat,Lon
    real, intent(out) :: interp_ktemp_space
    real, intent(out) :: interp_zonal_space
    real, intent(out) :: interp_merid_space
    real, intent(out) :: interp_density_space
    real, intent(out) :: interp_o_space
    real, intent(out) :: interp_o2_space
    real, intent(out) :: interp_nitrogen_space
    real, intent(out) :: interp_co2_space
    real, intent(out) :: interp_no_space
    real, intent(out) :: interp_n_space
    real, intent(out) :: interp_vert_space
 
    call interpolate_space(interp_ktemp,interp_ktemp_space)
    call interpolate_space(interp_zonal,interp_zonal_space)
    call interpolate_space(interp_merid,interp_merid_space)
    call interpolate_space(interp_density,interp_density_space)
    call interpolate_space(interp_o,interp_o_space)
    call interpolate_space(interp_o2,interp_o2_space)
    call interpolate_space(interp_nitrogen,interp_nitrogen_space)
    call interpolate_space(interp_co2,interp_co2_space)
    call interpolate_space(interp_no,interp_no_space)
    call interpolate_space(interp_atomic_n,interp_n_space)
    call interpolate_space(interp_vert,interp_vert_space)
   
    
    if (interp_ktemp_space<0) then
        write(*,*) 'negative temp in space'
    endif

    if (interp_density_space<0) then
        write(*,*) 'negative dens in space'
    endif

    if (interp_o_space<0) then
        write(*,*) 'negative o in space'
    endif

    if (interp_o2_space<0) then
        write(*,*) 'negative o2 in space'
    endif

    if (interp_nitrogen_space<0) then
        write(*,*) 'negative n2 in space'
    endif

    if (interp_co2_space<0) then
        write(*,*) 'negative co2 in space'
    endif

    if (interp_no_space<0) then
        write(*,*) 'negative no in space'
    endif

    if (interp_n_space<0) then
        write(*,*) 'negative n in space'
    endif

    contains
        subroutine interpolate_space(arr,arr_out)
            integer ::i,j,k=0
            integer :: i_l, i_u, j_l, j_u, k_l, k_u=0
            real::xd,yd,zd,c000,c100,c001,c101,c010,c110,c011,c111,c00,c01,c10,c11,c0,c1,c=0
            integer ::nLons,nLats,nAlts,nHours
            real, intent(in) ::arr(:,:,:)
            real, intent(out) ::arr_out
            integer,dimension(3)::shape_arr
            real :: dlon_waccm, dlat_waccm, dalt_waccm, mod_Lon

            shape_arr=shape(arr)
            nLons=shape_arr(1)
            nLats=shape_arr(2)
            nAlts=shape_arr(3)
            mod_Lon=modulo(Lon,360.0)! modulo of lon

            dlon_waccm=lon_waccm(2)-lon_waccm(1)
            dlat_waccm=lat_waccm(2)-lat_waccm(1)
            dalt_waccm=alt_waccm(2)-alt_waccm(1)
 
            do i=1,nLons
                if (lon_waccm(i)>mod_Lon) then
                    exit
                end if
            end do

            do j=1,nLats
                if (lat_waccm(j)>Lat) then
                    exit
                end if
            end do

            do k=1,nAlts
                if (alt_waccm(k)>Alt) then
                    exit
                end if
            end do
            
            ! maybe found is not required and lat< and lat>
            ! Ghost cell boundary conditions
            ! Circular interpolation for longitude
            if (i==1 .or. i==nLons+1) then
                 i_l=nLons
                 i_u=1
            else 
                i_l=i-1
                i_u=i
            end if
            
            if (j==1) then ! if Lat < -90.0
                j_l=1
                j_u=2
                yd=0
            else if (j==nLats+1) then ! if Lat > 90.0
                j_l=nLats-1
                j_u=nLats
                yd=1
            else
                j_l=j-1
                j_u=j
                yd=(Lat-lat_waccm(j_l))/dlat_waccm
            end if
            ! Coz this situation might arise only for lons and lats
            ! for alts not important because this code is being run
            ! for lower boundary conditions and gitm alts lie in 
            ! range 95, 97.5, 100 km
            k_l=k-1
            k_u=k
           
            xd=(mod_Lon-lon_waccm(i_l))/dlon_waccm
            zd=(Alt-alt_waccm(k_l))/dalt_waccm
            
            if (abs(xd)>1) then
                xd=((360.0+mod_Lon)-lon_waccm(i_l))/dlon_waccm
            end if

            !write(*,*) 'indices', i ,j,k
            !write(*,*) 'Grid Lon',Lon,lon_waccm(i),lon_waccm(i-1)
            !write(*,*) 'Grid Lat',Lat,lat_waccm(j), lat_waccm(j-1)
            !write(*,*) 'Grid Alt',Alt,alt_waccm(k), alt_waccm(k-1)
            !write(*,*) 'Coefs', xd,yd,zd
	    !if(minval(arr)<0.0) then
		!write(*,*) 'Starts arr negative'
	    !endif

            c000=arr(i_l,j_l,k_l)
            c100=arr(i_u,j_l,k_l)
            c001=arr(i_l,j_l,k_u)
            c101=arr(i_u,j_l,k_u)
            c010=arr(i_l,j_u,k_l)
            c110=arr(i_u,j_u,k_l)
            c011=arr(i_l,j_u,k_u)
            c111=arr(i_u,j_u,k_u)
 
            c00=(c000*(1-xd))+(c100*xd)
            c01=(c001*(1-xd))+(c101*xd)
            c10=(c010*(1-xd))+(c110*xd)
            c11=(c011*(1-xd))+(c111*xd)

            c0=(c00*(1-yd))+(c10*yd)
            c1=(c01*(1-yd))+(c11*yd)
                        
            c=(c0*(1-zd))+(c1*zd)
            arr_out=c

            !if (arr_out<=0.0) then
	    !write(*,*) 'Ends',i,j,k
	    !write(*,*) 'GITM',Lon,Lat,Alt
	    !write(*,*) 'Waccm',lon_waccm(i),lat_waccm(j),alt_waccm(k)
	    !write(*,*) 'Waccm previous',lon_waccm(i-1),lat_waccm(j-1),alt_waccm(k-1)
	    !write(*,*) 'Diffs', xd,yd,zd
            !write(*,*) 'c',c000,c100,c001,c101,c010,c110,c011,c111
	    !endif              
        end subroutine interpolate_space
end subroutine interp_3Dspace

subroutine interp_time
    use ModLowerBC 
    use ModTime, ONLY: CurrentTime
    use ModInputs
    implicit none
    call interpolate(ktemp,interp_ktemp)
    call interpolate(zonal,interp_zonal)
    call interpolate(merid, interp_merid)
    call interpolate(neutral_density,interp_density)
    call interpolate(o, interp_o)
    call interpolate(o2, interp_o2)
    call interpolate(nitrogen, interp_nitrogen)
    call interpolate(co2, interp_co2)
    call interpolate(no, interp_no)
    call interpolate(atomic_n, interp_atomic_n)
    call interpolate(vert,interp_vert)
    
    !write(*,*) 'Minval time : ', minval(interp_ktemp),minval(interp_atomic_n)
    if (minval(interp_ktemp)<0) then
        write(*,*) 'negative time temp'
    endif
    if (minval(atomic_n)<0)then
        write(*,*) 'First Minimum in time', minval(atomic_n)    
    endif
    
    if (minval(interp_atomic_n)<0) then
    write(*,*) 'Second Minimum in time', minval(interp_atomic_n)    
    endif

    
    !write(*,*) 'Hour_array', hour

    contains 
        subroutine interpolate(arr, arr_out)
    
            integer :: i, j, k, closest_index=0
            integer :: l=0
            real :: y1, y2, xval=0
            integer:: nLons, nLats, nAlts, nHours
            integer,dimension(4):: shape_arr
            real,intent(in):: arr(:,:,:,:)
            real,intent(out):: arr_out(:,:,:)
            double precision :: hour_UT, CurrentTime_UT
            double precision :: modulus
            integer :: real_hours
            !integer, dimension(7) :: iStartTime

            shape_arr=shape(arr)
            nLons=shape_arr(1)
            nLats=shape_arr(2)
            nAlts=shape_arr(3)
            nHours=shape_arr(4)
            
            !real_hours=nHours*nBCFiles
            real_hours=nHours

            !write(*,*) 'real',real_hours,nBCFiles,nLons,nLats,nAlts,nHours,CurrentTime
            !write(*,*) 'Current', CurrentTime

            !call time_real_to_int(CurrentTime,itime)
            !write(*,*) 'Beg', hour(1) 
            !write(*,*) 'shape', shape(hour)
            do l=1,real_hours
                if (hour(l) > CurrentTime) then
                    exit
                end if
            end do
             
            xval=(CurrentTime-hour(l-1))/(hour(l)-hour(l-1))
            !write(*,*) 'Interpolate Time',hour(l),CurrentTime, hour(l-1),&
            !            real_hours, xval
            do i=1,nLons
                do j=1,nLats
                    do k=1,nAlts
                        y1=arr(i,j,k,l-1)
                        y2=arr(i,j,k,l)
                        arr_out(i,j,k)=y1*(1-xval)+y2*(xval) 
                        !write(*,*) 'arr_out time : ', arr_out(i,j,k)
                    end do   
                end do
            end do
        end subroutine interpolate
end subroutine interp_time

subroutine init_waccm
    use ModGITM
    use ModVertical,only:Altitude_G
    use ModPlanet
    use ModConstants, only : Boltzmanns_Constant,pi
    use ModInputs

    implicit none

    integer :: iAlt,iBlock,iLon,iLat=0
    real :: msis_temp,g,geo_lat,geo_lon,dz,msis_temp_lower,Alt0,Alt,iSpecies
    real :: interp_ktemp_space, interp_zonal_space,&
        interp_merid_space,interp_density_space,interp_o_space,&
        interp_o2_space,interp_nitrogen_space,interp_co2_space,&
        interp_no_space,interp_n_space,interp_vert_space
    real :: scale_height_O,scale_height_O2,scale_height_N2,&
        scale_height_N, scale_height_NO
    
    do iBlock=1,nBlocks
        do iLon=-1,nLons+2
            do iLat=-1,nLats+2
                do iAlt=-1,0
                    
                    Alt0=Altitude_GB(iLon,iLat,iAlt,iBlock)/1000.0
                    geo_lat=Latitude(iLat,iBlock)*180.0/pi
                    geo_lon=mod(Longitude(iLon,iBlock)*180/pi+&
                        360.0,360.0)

                    call interp_time()
                    !write(*,*) 'Min val in O : ', minval(interp_o)
                    !write(*,*) 'Min val in O2 : ', minval(interp_o2)
                    !write(*,*) 'Time Done'
                    call interp_3Dspace(Alt0,&
                        geo_lat,geo_lon,&
                        interp_ktemp_space, interp_zonal_space,&
                        interp_merid_space,interp_density_space,&
                        interp_o_space, interp_o2_space,&
                        interp_nitrogen_space,interp_co2_space,&
                        interp_no_space,interp_n_space,interp_vert_space)
                    !write(*,*) 'Space Done'
                    if (.not. LowBMSISonlyO) then
                        NDensityS(iLon,iLat,iAlt,iO_3P_,iBlock)=&
                            interp_o_space                    
                        if (interp_o_space<=0.0) then
			    write(*,*) 'Here it is 1'                    
                        endif
                    endif

		    NDensityS(iLon,iLat,iAlt,iO2_,iBlock)=&
                        interp_o2_space
                    if (interp_o2_space<=0.0) then
			write(*,*) 'Here it is 2'                    
                    endif

		    NDensityS(iLon,iLat,iAlt,iN2_,iBlock)=&
                        interp_nitrogen_space                    
                    if (interp_nitrogen_space<=0.0) then
			write(*,*) 'Here it is 3'                    
                    endif


                    !NDensityS(iLon,iLat,iAlt,iN_4S_,iBlock)=&
                    !    interp_n_space
		    !if (interp_n_space<=0.0) then
			!write(*,*) 'Here it is 4'                    
                    !endif
		    
		    NDensityS(iLon,iLat,iAlt,iNO_,iBlock)=&
                        interp_no_space
		    if (interp_no_space<=0.0) then
			write(*,*) 'Here it is 5'                    
                    endif
                    
                    !write(*,*) 'Hello 2 :',interp_o_space

                    !write(*,*) 'Temp : ', interp_ktemp_space, &
                    !    TempUnit(iLon,iLat,iAlt),iLon,iLat,iAlt,iBlock
                    Temperature(iLon,iLat,iAlt,iBlock)=interp_ktemp_space/&
                        TempUnit(iLon,iLat,iAlt)
                    
                    Velocity(iLon,iLat,iAlt,iEast_,iBlock)=interp_zonal_space
                    Velocity(iLon,iLat,iAlt,iNorth_,iBlock)=interp_merid_space
                    Velocity(iLon,iLat,iAlt,iUp_,iBlock)=interp_vert_space

                enddo
                do iAlt=1, nAlts+2
                    dz=(Altitude_GB(iLon,iLat,iAlt,iBlock)-&
                        Altitude_GB(iLon,iLat,iAlt-1,iBlock))/1000.0
                 
                    msis_temp=Temperature(iLon,iLat,iAlt,iBlock)*&
                        TempUnit(iLon,iLat,iAlt)
                    msis_temp_lower=Temperature(iLon,iLat,iAlt-1,iBlock)*&
                        TempUnit(iLon,iLat,iAlt)

                    g=Gravity_GB(iLon,iLat,iAlt,iBlock)
                    scale_height_O=-(Boltzmanns_Constant*msis_temp/&
                        (g*Mass(1)))/1000.0
                    scale_height_O2=-(Boltzmanns_Constant*msis_temp/&
                        (g*Mass(2)))/1000.0
                    scale_height_N2=-(Boltzmanns_Constant*msis_temp/&
                        (g*Mass(3)))/1000.0
                    scale_height_N=-(Boltzmanns_Constant*msis_temp/&
                        (g*Mass(4)))/1000.0
                    scale_height_NO=-(Boltzmanns_Constant*msis_temp/&
                        (g*Mass(5)))/1000.0
                
                    Alt=Altitude_GB(iLon,iLat,iAlt,iBlock)/1000.0
                    if (.not. LowBMSISonlyO) then
                        NDensityS(iLon,iLat,iAlt,iO_3P_,iBlock)=&
                            NDensityS(iLon,iLat,iAlt-1,iO_3P_,iBlock)&
                            *exp(-dz/scale_height_O)*(msis_temp_lower/msis_temp)
                        if(minval(NDensityS(iLon,iLat,:,iO_3P_,iBlock))<=0.0) then
                            write(*,*) "negative density found in 1 !"
                        endif
                    endif

                    NDensityS(iLon,iLat,iAlt,iO2_,iBlock)=&
                    NDensityS(iLon,iLat,iAlt-1,iO2_,iBlock)&
                    *exp(-dz/scale_height_O2)*(msis_temp_lower/msis_temp)

                    if(minval(NDensityS(iLon,iLat,:,iO2_,iBlock))<=0.0) then
                            write(*,*) "negative density found in 2 !"
                    endif

                    NDensityS(iLon,iLat,iAlt,iN2_,iBlock)=&
                    NDensityS(iLon,iLat,iAlt-1,iN2_,iBlock)&
                    *exp(-dz/scale_height_N2)*(msis_temp_lower/msis_temp)

                    if(minval(NDensityS(iLon,iLat,:,iN2_,iBlock))<=0.0) then
                            write(*,*) "negative density found in 3 !"
                    endif

                    !NDensityS(iLon,iLat,iAlt,iN_4S_,iBlock)=&
                    !NDensityS(iLon,iLat,iAlt-1,iN_4S_,iBlock)&
                    !*exp(-dz/scale_height_N)*(msis_temp_lower/msis_temp)

                    !if(minval(NDensityS(iLon,iLat,:,iN_4S_,iBlock))<=0.0) then
                    !	write(*,*) 'Minimum Density', minval(NDensityS(iLon,iLat,:,iN_4S_,iBlock))
                    !    write(*,*) "negative density found in 4 !"
                    
                    !endif

                    NDensityS(iLon,iLat,iAlt,iNO_,iBlock)=&
                    NDensityS(iLon,iLat,iAlt-1,iNO_,iBlock)&
                    *exp(-dz/scale_height_NO)*(msis_temp_lower/msis_temp)
                    
                    if(minval(NDensityS(iLon,iLat,:,iNO_,iBlock))<=0.0) then
                            write(*,*) "negative density found in 5 !"
                    endif

                    

                    MeanMajorMass(iLon,iLat,iAlt)=0
                    do iSpecies=1, nSpecies
                        MeanMajorMass(iLon,iLat,iAlt)=&
                        MeanMajorMass(iLon,iLat,iAlt)+&
                        Mass(iSpecies)*NDensityS(iLon,iLat,&
                        iAlt,iSpecies,iBlock)/&
                        sum(NDensityS(iLon,iLat,iAlt,&
                        1:nSpecies,iBlock))
                    enddo


                    LogNS(iLon,iLat,iAlt,1:nSpecies,iBlock)=&
                    log(NDensityS(iLon,iLat,iAlt,1:nSpecies,iBlock))
                     
                enddo
            enddo
        enddo

    !write(*,*) 'hello4'    
    Rho(:,:,:,iBlock)=0.0
    NDensity(:,:,:,iBlock)=0.0

    do iSpecies=1,nSpecies
        NDensity(:,:,:,iBlock)=NDensity(:,:,:,iBlock)+&
            NDensityS(:,:,:,iSpecies,iBlock)
        Rho(:,:,:,iBlock)=Rho(:,:,:,iBlock)+&
            Mass(iSpecies)*NDensityS(:,:,:,iSpecies,iBlock)
    enddo   

    enddo
    
end subroutine init_waccm






