  
!===============================================================================
 
 subroutine OA_ncreadstat(nc_fname_fld, varapr_fname, covsca_fname,     &
                          bathy_fname,                                  &
			  PARAM, fact, iarea, nx_ana, ny_ana, nb_level, &
			  nx_bathy, ny_bathy,                           &
                          nc_status, lat_ana, lon_ana, ana_3Dvartot,    &
                          covar_ms, bath_ana, msk_area,                 &
                          lat_bat, lon_bat, covar_all, bath_all)

!===============================================================================
!
!  Objet: Reads field grid and covariances in nc file prepared by pre-OA
!  -----    
!
!  Version: 4.0 f95 - mpi
!  -------
! history:
! V4.0	: 22/05/2006: 	F. Gaillard	- F95-MPI
! V4.01	: 15/05/2007: 	F. Gaillard	- F95-
!
! Calls: 
!	netcdf library
!		
!===============================================================================
!

use netcdf

implicit none


!  INPUTS:
!  ------
character*(*) :: nc_fname_fld, varapr_fname, covsca_fname, bathy_fname
character*4 PARAM
real*8    :: fact
integer*4 :: iarea, nx_ana, ny_ana, nb_level, nx_bathy, ny_bathy 

!  OUTPUTS:
!  ------
integer*4 :: nc_status      ! status index after reading nc_file 
real*8    :: lat_ana(*), lon_ana(*)
real*8    :: ana_3Dvartot(nx_ana, ny_ana, nb_level)	! Total variance 
real*8    :: covar_ms(nx_ana,ny_ana,2)	! meso  scale covariances (X and Y)
real*8    :: bath_ana(nx_ana,ny_ana)	
integer*4 :: msk_area(nx_ana,ny_ana)
real*8    :: lat_bat(*), lon_bat(*)
real*8    :: covar_all(nx_bathy, ny_bathy,2), &  ! bathy and covariances over global ocean
             bath_all(nx_bathy, ny_bathy)

!  LOCALS:
!  ------
integer*4    :: lu_nc, dim_id, id_var, i, j, ii, jj, ij, i0, j0 , &
                area_typ, fill_val, ln1
real*4       :: xoff, scal_fact
character*8  :: par_std_name

!  DYNAMICALLY ALLOCATED ARRAYS

real*4 ,   allocatable :: tab1Dx(:), tab1Dy(:), cov_ms(:,:)
integer*4, allocatable :: ii0(:), jj0(:)
integer*2, allocatable :: tab_i2(:,:), std_3D(:,:,:), tab_out(:,:,:)


!  --------------------------------
!  Program starts
!  --------------------------------

par_std_name = PARAM//'_STD'

! ====================================
!    Coordinates of analysis points:
! ====================================

 ! Opens nc file 
 lu_nc     = 30
 ln1 = index(nc_fname_fld, '.', back=.true.) + 2 
 nc_status = nf90_open(nc_fname_fld(1:ln1), nf90_nowrite, lu_nc)
 if (nc_status .ne. 0) then
    nc_status = 100+nc_status
    return
 end if

 allocate(tab1Dy(ny_ana), tab1Dx(nx_ana))

 nc_status = nf90_inq_varid(lu_nc, 'longitude', id_var)
 nc_status = nf90_get_var(lu_nc, id_var, tab1Dx)
 nc_status = nf90_inq_varid(lu_nc, 'latitude', id_var)
 nc_status = nf90_get_var(lu_nc, id_var, tab1Dy)
 nc_status = nf90_close (lu_nc)

 lat_ana(1:ny_ana) = dble(tab1Dy)
 lon_ana(1:nx_ana) = dble(tab1Dx)



 deallocate(tab1Dx, tab1Dy)



!  =================
!     Bathymetry :
!  =================
 ! Opens nc file 
 lu_nc     = 30
 ln1 = index(bathy_fname, '.', back=.true.) + 2 
 nc_status = nf90_open(bathy_fname(1:ln1), nf90_nowrite, lu_nc)
 if (nc_status .ne. 0) then
    nc_status = 200+nc_status
    return
 end if

!  Reads coordinates
!  ------------------
 allocate(tab1Dy(ny_bathy), tab1Dx(nx_bathy))
 nc_status = nf90_inq_varid(lu_nc, 'longitude', id_var)
 nc_status = nf90_get_var(lu_nc, id_var, tab1Dx)

 nc_status = nf90_inq_varid(lu_nc, 'latitude', id_var)
 nc_status = nf90_get_var(lu_nc, id_var, tab1Dy)

 lon_bat(1:nx_bathy) = dble(tab1Dx)
 lat_bat(1:ny_bathy) = dble(tab1Dy)
 deallocate(tab1Dx, tab1Dy)


!  Identifies indices of area in global table
!  ------------------------------------------
 allocate(ii0(ny_bathy), jj0(nx_bathy))
 ii0 = minloc(abs(lat_bat(1:ny_bathy) - lat_ana(1))) 
 jj0 = minloc(abs(lon_bat(1:nx_bathy) - lon_ana(1))) 
 i0 = ii0(1)
 j0 = jj0(1)
 deallocate(ii0, jj0)
 
!  Reads full bathymetry
!  ---------------------
 allocate(tab_i2(nx_bathy,ny_bathy))
 nc_status = nf90_inq_varid(lu_nc, 'bathymetry', id_var)
 nc_status = nf90_get_var(lu_nc, id_var, tab_i2)
 bath_all  = dble(tab_i2)
 deallocate(tab_i2)

!  Extracts local bathymetry
!  -------------------------
 bath_ana = bath_all(j0:j0-1+nx_ana, i0:i0-1+ny_ana)
 



!  Reads local mask
!  ---------------------
 allocate(tab_i2(nx_ana,ny_ana))					
 nc_status = nf90_inq_varid(lu_nc, 'basin_area', id_var)
 nc_status = nf90_get_var(lu_nc, id_var, tab_i2, start = (/j0, i0/), &
                                                count = (/nx_ana,ny_ana/))
 msk_area = tab_i2
 deallocate(tab_i2)

 nc_status = nf90_close (lu_nc)

 
!  ========================
!    Covariance scales:
!  ========================

! Opens nc file 
 lu_nc     = 30    
 ln1 = index(covsca_fname, '.') + 2 
 nc_status = nf90_open(covsca_fname(1:ln1), nf90_nowrite, lu_nc)
 if (nc_status .ne. 0) then
    nc_status = 300+nc_status
    return
 end if

 allocate(cov_ms(nx_bathy,ny_bathy))

 nc_status = nf90_inq_varid(lu_nc, 'COV_SCALE', id_var)
 if (nc_status .eq. 0) then !old format
     nc_status = nf90_get_var(lu_nc, id_var, cov_ms)
     covar_all(:,:,1) = dble(cov_ms)
     covar_ms(:,:,1) = covar_all(j0:j0-1+nx_ana, i0:i0-1+ny_ana,1)
     covar_all(:,:,2) = covar_all(:,:,1)
     covar_ms(:,:,2) =  covar_ms(:,:,1)

 else  !new format   
     nc_status = nf90_inq_varid(lu_nc, 'COV_SCALE_X', id_var)
     nc_status = nf90_get_var(lu_nc, id_var, cov_ms)
     covar_all(:,:,1) = dble(cov_ms)
     covar_ms(:,:,1) = covar_all(j0:j0-1+nx_ana, i0:i0-1+ny_ana,1)
 
     nc_status = nf90_inq_varid(lu_nc, 'COV_SCALE_Y', id_var)
     nc_status = nf90_get_var(lu_nc, id_var, cov_ms)
     covar_all(:,:,2) = dble(cov_ms)
     covar_ms(:,:,2) = covar_all(j0:j0-1+nx_ana, i0:i0-1+ny_ana,2)
 
 end if

 nc_status = nf90_close (lu_nc)

 deallocate (cov_ms)
 
 
 
!  ======================
!     Variances:
!  ====================== 

 nc_status = nf90_open(varapr_fname, nf90_nowrite, lu_nc)
 if (nc_status .ne. 0) then
    nc_status = 400+nc_status
    return
 end if

 allocate(std_3D(nx_ana, ny_ana, nb_level), tab_out(nx_ana, ny_ana, nb_level))

 
 nc_status = nf90_inq_varid(lu_nc, par_std_name, id_var)
 nc_status = nf90_get_var(lu_nc, id_var, std_3D, start = (/j0, i0, 1/), &
                                    count = (/nx_ana,ny_ana,nb_level/))
 nc_status = nf90_get_att(lu_nc, id_var, 'add_offset', xoff)	
 nc_status = nf90_get_att(lu_nc, id_var, 'scale_factor', scal_fact)
 nc_status = nf90_get_att(lu_nc, id_var, '_FillValue', fill_val)		    
 nc_status = nf90_close (lu_nc)
 tab_out(1:nx_ana,1:ny_ana,1:nb_level) = 1
 where(std_3D(1:nx_ana,1:ny_ana,1:nb_level) .EQ. fill_val) &
      tab_out(1:nx_ana,1:ny_ana,1:nb_level) = 0
 
 ana_3Dvartot = xoff + dble(std_3D)*scal_fact
 ana_3Dvartot = fact*ana_3Dvartot*ana_3Dvartot*tab_out
 
! write(*,*) 'xoff, scal_fact, fill_val'
! write(*,*) xoff, scal_fact, fill_val
! write (*,*) 'ana_3Dvartot(1,1,:)'
! write (*,*) ana_3Dvartot(1,1,:)
! write (*,*) 'ana_3Dvartot(1,:,1)'
! write (*,*) ana_3Dvartot(1,:,1)
! write (*,*) 'ana_3Dvartot(:,1,1)'
! write (*,*) ana_3Dvartot(:,1,1)
 
 deallocate(std_3D, tab_out)
 
 ! adjust mask with area where coviance is defined and bathy > 4m
 where(covar_ms(1:nx_ana,1:ny_ana,1) .LE. 0)  msk_area(1:nx_ana,1:ny_ana) = -1
 where(covar_ms(1:nx_ana,1:ny_ana,2) .LE. 0)  msk_area(1:nx_ana,1:ny_ana) = -1
 where(ana_3Dvartot(1:nx_ana,1:ny_ana,1) .LE. 0)  msk_area(1:nx_ana,1:ny_ana) = -1
 where(bath_ana(1:nx_ana,1:ny_ana) .LT. 4)    msk_area(1:nx_ana,1:ny_ana) = -1


!write(*,*) lat_ana(1:ny_ana)
!write(*,*) lon_ana(1:nx_ana)

 !write (*,*) 'bath_ana(1,1:ny_ana)'
 !write (*,*) bath_ana(1,1:ny_ana)
 !write (*,*) 'bath_ana(1:nx_ana,1)'
 !write (*,*) bath_ana(1:nx_ana,1)

 !write (*,*) 'ana_3Dvartot (1,1:ny_ana,1)'
 !write (*,*) ana_3Dvartot(1,1:ny_ana,1)
 !write (*,*) 'ana_3Dvartot (1:nx_ana,1,1)'
 !write (*,*) ana_3Dvartot(1:nx_ana,1,1)

 !write(*,*) 'covar_ms(1,1:ny_ana,1'
 !write (*,*) covar_ms(1,1:ny_ana,1)
 !write(*,*) 'covar_ms(1:nx_ana,1,1'
 !write (*,*) covar_ms(1:nx_ana,1,1)

 !write(*,*) 'covar_ms(1,1:ny_ana,2'
 !write (*,*) covar_ms(1,1:ny_ana,2)
 !write(*,*) 'covar_ms(1:nx_ana,1,2'
 !write (*,*) covar_ms(1:nx_ana,1,2)

 return
 end
