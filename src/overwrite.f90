
!---------------------------------------------------------------------------
! Routines for making GITM use another model of the ionosphere
!---------------------------------------------------------------------------

subroutine overwrite_ionosphere

  use ModInputs
  use ModSizeGitm
  use ModGITM, only : iProc
  
  implicit none

  integer :: iBlock

  if (DoOverwriteWithIRI) then
     if (iProc == 0) write(*,*) 'Overwriting Ionosphere with IRI!'
     ! This doesn't include the velocities
     call init_iri
  else
     if (iProc == 0) write(*,*) 'Overwriting Ionosphere with SAMI-3!'
     do iBlock = 1, nBlocks
        call ionosphere_overwrite_sami(iBlock)
     enddo
  endif
  
end subroutine overwrite_ionosphere

!---------------------------------------------------------------------------
! Routine for overwriting GITM data with SAMI data
!---------------------------------------------------------------------------

subroutine ionosphere_overwrite_sami(iBlock)

  use ModReadSami3D

  use ModSizeGitm
  use ModPlanet, only : nSpecies, nIons
  use ModInputs
  use ModReadSami3d
  use ModGITM
  use ModTime, only : CurrentTime
  
  implicit none

  save
  
  integer, intent(in) :: iBlock
  character (len=nSamiVarCharLength), allocatable :: SamiVars(:)
  real, allocatable :: GitmLons(:),GitmLats(:),GitmAlts(:)
  real, allocatable :: SamiFileData(:,:)
  logical :: IsFirstTime = .true.
  
  integer :: nPoints, nVars, iErr
  integer :: iPoint, iLon, iLat, iAlt

  call report("ionosphere_overwrite_sami",2)

  if (IsFirstTime) then

     iErr = 0

     call SamiReadInputFile(SamiInFile)
     call SamiGetnVars(nVars)

     nPoints = (nLons+4) * (nLats+4) * (nAlts)
     call SamiSetnPointsToGet(nPoints)
     
     allocate(SamiVars(nVars))
     allocate(SamiFileData(nPoints,nVars))

     call SamiGetVars(SamiVars)

     allocate(GitmLons(nPoints))
     allocate(GitmLats(nPoints))
     allocate(GitmAlts(nPoints))

     iPoint = 1

     do iLon = -1,nLons+2
        do iLat = -1, nLats+2
           do iAlt = 1, nAlts
              GitmLons(iPoint) = longitude(iLon,iBlock)
              GitmLats(iPoint) = latitude(iLat,iBlock)
              GitmAlts(iPoint) = Altitude_GB(iLon,iLat,iAlt,iBlock)
              iPoint = iPoint + 1
           enddo
        enddo
     enddo

     if (iPoint > 1) then 

        ! Convert to degrees and km
        GitmLons = GitmLons * 360.0 / twopi
        GitmLats = GitmLats * 360.0 / twopi
        GitmAlts = GitmAlts / 1000.0
     
        call SamiSetGrid(GitmLons,GitmLats,GitmAlts)

     endif
     
     IsFirstTime = .false.
     
  endif

  if (CorotationAdded) then 
     ! This then implies that the grid is in local time, so we need to 
     ! feed SAMI the local time (still 0-360, though).

     iPoint = 1

     do iLon = -1,nLons+2
        do iLat = -1, nLats+2
           do iAlt = 1, nAlts
              GitmLons(iPoint) = LocalTime(iLon)*15.0
              if (iAlt == 1 .and. iLat == 1) write(*,*) iLon, LocalTime(iLon), GitmLons(iPoint) 
              iPoint = iPoint + 1
           enddo
        enddo
     enddo

     if (iPoint > 1) call SamiSetGrid(GitmLons,GitmLats,GitmAlts)

  endif

  call SamiUpdateTime(CurrentTime, iErr)
  if (iErr == 0) then 
     call SamiGetData(SamiFileData)
  endif

  iPoint = 1

  do iLon = -1,nLons+2
     do iLat = -1, nLats+2
        do iAlt = 1, nAlts
           !write(*,*) iLon,iLat,iAlt,iPoint
           IDensityS(iLon,iLat,iAlt,iO_4SP_,iBlock) = SamiFileData(iPoint, iSami_Op_)*1e6
           IDensityS(iLon,iLat,iAlt,iO2P_,iBlock) = SamiFileData(iPoint, iSami_O2p_)*1e6
           IDensityS(iLon,iLat,iAlt,iNOP_,iBlock) = SamiFileData(iPoint, iSami_NOp_)*1e6
           IDensityS(iLon,iLat,iAlt,iHP_,iBlock) = SamiFileData(iPoint, iSami_Hp_)*1e6
           IDensityS(iLon,iLat,iAlt,iHeP_,iBlock) = SamiFileData(iPoint, iSami_Hep_)*1e6

           IDensityS(iLon,iLat,iAlt,ie_,iBlock) = &
                sum(IDensityS(iLon,iLat,iAlt,1:nIons-1,iBlock))

           eTemperature(iLon,iLat,iAlt,iBlock) = SamiFileData(iPoint, iSami_Te_)
           iTemperature(iLon,iLat,iAlt,iBlock) = SamiFileData(iPoint, iSami_Ti_)

           iPoint = iPoint + 1
        enddo
     enddo
  enddo
  
  
end subroutine ionosphere_overwrite_sami

subroutine overwrite_thermosphere_winds 
     !need to overwrite thermospheric winds
     use ModGITM !,only : Velocity,longitude,latitude,Altitude_GB
     use ModTime ,only : tSimulation
     use ModInputs
     use ModSizeGitm
     use ModConstants, only: pi
     real :: interp_ktemp_space,interp_zonal_space,interp_merid_space,&
      interp_density_space,interp_o_space,interp_o2_space,&
      interp_nitrogen_space,interp_co2_space,interp_no_space,&
      interp_n_space,interp_vert_space
     !real, dimension(1:nAlts):: vert_prof
     real :: vert_prof,alpha,zlb,zmax,temp_north, temp_east, temp_vert

     !TimeConstant=5*60
     !Timefactor=exp(-tSimulation/TimeConstant)
     G=1/RelaxationTime
     alpha=G*Dt
     alpha=1
     !write(*,*) ,'Relaxation : ', RelaxationTime, G, Dt,alpha
     zlb=95
     zmax=140
     do iBlock = 1, nBlocks
        do iLon = -1,nLons+2
            do iLat = -1, nLats+2
                do iAlt = 1, nAlts
                    gitm_lat=Latitude(iLat,iBlock)*180.0/pi
                    gitm_lon=mod(Longitude(iLon,iBlock)*180/pi+&
                        360.0,360.0)
                    gitm_alt= Altitude_GB(iLon,iLat,iAlt,iBlock)/1000.0
                    frac=(gitm_alt-zlb)/(zmax-zlb)
                    vert_prof=(cos((pi/2)*frac))**1.2
                    if (gitm_alt<=zmax) then
                        call interp_3Dspace(gitm_alt,gitm_lat,gitm_lon,interp_ktemp_space,interp_zonal_space,&
                            interp_merid_space,interp_density_space,interp_o_space,interp_o2_space,&
                            interp_nitrogen_space,interp_co2_space,interp_no_space,&
                            interp_n_space,interp_vert_space)
                            
                            !write(*,*) 'Indices : ', iBlock, iLon, iLat, iAlt
                            !write(*,*)  'GITM alt:',gitm_alt, gitm_lat, gitm_lon

                            !write(*,*) 'Alt : ', gitm_alt,  'Before Merid : ',Velocity(iLon,iLat,iAlt,iNorth_,iBlock),'Before Zonal : ', Velocity(iLon,iLat,iAlt,iEast_,iBlock),&
                            !'vert_prof : ', vert_prof,'interp merid : ', interp_merid_space, 'interp zonal : ', interp_zonal_space
                            
                            temp_north =((1-alpha*vert_prof)* &
                                Velocity(iLon,iLat,iAlt,iNorth_,iBlock))+((alpha*vert_prof)*interp_merid_space)
                            temp_east=((1-alpha*vert_prof)* &
                                Velocity(iLon,iLat,iAlt,iEast_,iBlock))+((alpha*vert_prof)*interp_zonal_space)
                            temp_vert=((1-alpha*vert_prof)* &
                                Velocity(iLon,iLat,iAlt,iUp_,iBlock))+((alpha*vert_prof)*interp_vert_space)

                            Velocity(iLon,iLat,iAlt,iNorth_,iBlock)=temp_north
                            Velocity(iLon,iLat,iAlt,iEast_,iBlock)=temp_east
                            Velocity(iLon,iLat,iAlt,iUp_,iBlock)=temp_vert

                            !write(*,*) 'Alt : ', gitm_alt, 'After Merid : ',temp_north, Velocity(iLon,iLat,iAlt,iNorth_,iBlock),' After Zonal : ', temp_east, Velocity(iLon,iLat,iAlt,iEast_,iBlock),&
                            !'vert_prof : ',vert_prof,'interp merid : ', interp_merid_space, 'interp zonal : ', interp_zonal_space

                            !Velocity(iLon,iLat,iAlt,iUp_,iBlock)=((1-alpha*vert_prof)* &
                            !    Velocity(iLon,iLat,iAlt,iNorth_,iBlock))+((alpha*vert_prof)*interp_merid_space)
                            !Temperature(iLon,iLat,iAlt,iBlock)=((1-alpha*vert_prof)* &
                            !    Temperature(iLon,iLat,iAlt,iNorth_,iBlock))+((alpha*vert_prof)*interp_ktemp_space)
                    endif
                enddo
            enddo
        enddo
      enddo
    end subroutine overwrite_thermosphere_winds
