module ModLowerBC
    implicit none
    real,allocatable::lat_waccm(:),lon_waccm(:),&
        alt_waccm(:),hour(:)
    real,allocatable:: hour_arr(:)
    real,allocatable::ktemp(:,:,:,:)
    real,allocatable::zonal(:,:,:,:)
    real,allocatable::merid(:,:,:,:)
    real,allocatable::neutral_density(:,:,:,:)
    real,allocatable::o(:,:,:,:)
    real,allocatable::o2(:,:,:,:)
    real,allocatable::nitrogen(:,:,:,:)
    real,allocatable::co2(:,:,:,:)
    real,allocatable::no(:,:,:,:)
    real,allocatable::atomic_n(:,:,:,:)
    real,allocatable::vert(:,:,:,:)

    real,allocatable::interp_ktemp(:,:,:)
    real,allocatable::interp_zonal(:,:,:)
    real,allocatable::interp_merid(:,:,:)
    real,allocatable::interp_density(:,:,:)
    real,allocatable::interp_o(:,:,:)
    real,allocatable::interp_o2(:,:,:)
    real,allocatable::interp_nitrogen(:,:,:)
    real,allocatable::interp_co2(:,:,:)
    real,allocatable::interp_no(:,:,:)
    real,allocatable::interp_atomic_n(:,:,:)
    real,allocatable::interp_vert(:,:,:)

end module ModLowerBC
