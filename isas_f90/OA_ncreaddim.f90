 

!===============================================================================
 
   subroutine OA_ncreaddim(nc_fname_fld, nc_fname_dat, bathy_fname, &
                           nc_status, nb_profile, nx_ana, ny_ana, nb_level, & 
                           jjul_est, nx_bathy, ny_bathy)

!===============================================================================
!
!  Objet: Reads dimensions in nc file prepared by pre-OA
!  -----    
!
!  Version: 4.0 f95 - mpi
!  -------
! history:
! V4.0	: 07/06/2006: 	F. Gaillard	- F95-MPI
! V4.01	: 15/05/2007: 	F. Gaillard	- F95-
!        switch to nf90
! V4.01	: 15/05/2007: 	F. Gaillard	- Reads also bathymetry dimensions
! Calls: 
!	netcdf library
!		
!===============================================================================
!

use netcdf

implicit none



!  INPUTS:
!  ------
character*(*) :: nc_fname_fld, &  ! Name of nc file holding the field
                 nc_fname_dat, &  ! Name of nc file holding the data
                 bathy_fname      ! Name of nc file holding the bathymetry

!  OUTPUTS:
!  ------
integer*4 :: nc_status,         & ! status index after reading nc_file 
             nb_profile,        & ! Number of data profiles (nb_profile<= i_data)
	     nx_ana, ny_ana,    & ! Number of analysis points in x(long) and y (lat)
             nb_level,	        & ! Number of data levels (nb_level <= i_level)
             nx_bathy, ny_bathy   ! dimensions of global fields

real*8 :: jjul_est


!  LOCALS:
!  ------
integer*4    :: lu_nc, id_var, dim_id, ln1
real*4       :: jjul4


!  --------------------------------
!  Program starts
!  --------------------------------

 lu_nc     = 30    

!write(*,*) nc_fname_fld
! write(*,*) nc_fname_dat
! write(*,*) bathy_fname

! Opens field  file 
!  ------------------
 ln1 = index(nc_fname_fld, '.', back=.true.) + 2 
 nc_status = nf90_open(nc_fname_fld(1:ln1), nf90_nowrite, lu_nc)
 if (nc_status .ne. 0) then
     nc_status = 100 + nc_status
     return
 endif

!  Dimensions
 nc_status = nf90_inq_dimid(lu_nc, 'depth', dim_id)
 nc_status = nf90_Inquire_Dimension(lu_nc, dim_id, len = nb_level)
 if (nc_status .ne. 0) return
! write(*,*) nb_level
 
 nc_status = nf90_inq_dimid(lu_nc, 'longitude', dim_id)
 nc_status = nf90_Inquire_Dimension(lu_nc, dim_id, len = nx_ana)
 if (nc_status .ne. 0) return

 nc_status = nf90_inq_dimid(lu_nc, 'latitude', dim_id) 
 nc_status = nf90_Inquire_Dimension(lu_nc, dim_id, len = ny_ana)
 if (nc_status .ne. 0) return
 
 nc_status = nf90_inq_varid(lu_nc, 'time', id_var)
 nc_status = nf90_get_var(lu_nc, id_var, jjul4)
 if (nc_status .ne. 0) return
 jjul_est = dble(jjul4)
 
 nc_status = nf90_close (lu_nc)
 
 
! Opens data file 
!  ------------------
ln1 = index(nc_fname_dat, '.') + 2 
nc_status = nf90_open(nc_fname_dat(1:ln1), nf90_nowrite, lu_nc)
if (nc_status .ne. 0) return


!  Dimensions
 nc_status = nf90_inq_dimid(lu_nc, 'N_PROF', dim_id)
 nc_status = nf90_Inquire_Dimension(lu_nc, dim_id, len = nb_profile)
 if (nc_status .ne. 0) return
 
 nc_status = nf90_close (lu_nc)

! Opens bathy file 
!  ------------------
 ln1 = index(bathy_fname, '.') + 2 
 nc_status = nf90_open(bathy_fname(1:ln1), nf90_nowrite, lu_nc)
 if (nc_status .ne. 0) then
     nc_status = 200 + nc_status
     return
 endif


!  Dimensions
 nc_status = nf90_inq_dimid(lu_nc, 'longitude', dim_id)
 nc_status = nf90_Inquire_Dimension(lu_nc, dim_id, len = nx_bathy)
 nc_status = nf90_inq_dimid(lu_nc, 'latitude', dim_id)
 nc_status = nf90_Inquire_Dimension(lu_nc, dim_id, len = ny_bathy)
 nc_status = nf90_close (lu_nc)

return
end
