 
!===============================================================================
 
  subroutine  OA_ncwrite_fld (nc_fname, PARAM, nb_level, ny_ana, nx_ana, &
                              ana_3Dfld, ana_3Dpctvar, ana_3Dvarapr, nc_status )
 
!===============================================================================
!
!  Objet: Writes analysis results in nc file prepared by pre-OA
!  -----    
!
!  Version: 4.0 f95 - mpi
!  -------
! history:
! V4.0	: 22/05/2006: 	F. Gaillard	- F95-MPI
! V5.0	: 29/04/2009: 	F. Gaillard	- F95
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
 character*(*) :: nc_fname 	! Name of nc file holding the data
 character*4   :: PARAM 		! Type of variable (TEMP or PSAL)
 integer*4     :: nb_level,  &   ! number of levels  
                  ny_ana, nx_ana ! number of y (lat) and x(lon) points in analysis
		                ! grid
 real*8   :: ana_3Dfld(nx_ana, ny_ana, nb_level),   & ! Analyzed field (AF)
             ana_3Dpctvar(nx_ana, ny_ana, nb_level),& ! Percent variance on AF 
	     ana_3Dvarapr(nx_ana, ny_ana, nb_level)  ! A priori variance on AF 
				
!  OUTPUTS:
!  ------
 integer*4 :: nc_status    ! status index after reading nc_file 


!  LOCALS:
!  ------
 integer*4    :: lu_nc, io_mode, id_var
 real*4       :: xoff, scal_fact, xval, xmax
 integer*2    :: fill_val, ii, imax
		
 character*20   :: var_name

!  DYNAMICALLY ALLOCATED ARRAYS

 integer*2, allocatable :: std_3D(:,:,:)
 integer*1, allocatable :: pct_3D(:,:,:)

!  --------------------------------
!  Program starts
!  --------------------------------

! Opens nc file 
 lu_nc     = 30    
 nc_status = nf90_open(nc_fname, nf90_write, lu_nc)
 if (nc_status .ne. 0) return

 
 ! later save also: fact, qc_max, cov_LS, var_weigh, conf_xyzt, covar_ls, covar_ms_t,

 allocate (std_3D(nx_ana,ny_ana, nb_level))
 allocate (pct_3D(nx_ana,ny_ana, nb_level))
 
 !! copy analized field in nc file
 var_name = PARAM//'_ANO' !! attention only anomaly at this point
 nc_status = nf90_inq_varid(lu_nc, var_name, id_var)
 if (nc_status .ne. 0) return
 nc_status = nf90_get_att(lu_nc, id_var, 'scale_factor', scal_fact)
 nc_status = nf90_get_att(lu_nc, id_var, '_FillValue', fill_val)
 imax = abs(fill_val) - 1
 xmax = imax*scal_fact
 std_3D = ana_3Dfld /scal_fact 
 where (abs(ana_3Dfld) >= xmax)  std_3D = fill_val
 nc_status = nf90_put_var(lu_nc, id_var, std_3D)
 if (nc_status .ne. 0) return
  
 !! copy pctvar field in nc file
 var_name = PARAM//'_ANO_PCTVAR'
 nc_status = nf90_inq_varid(lu_nc, var_name, id_var)
 if (nc_status .ne. 0) return
 nc_status = nf90_get_att(lu_nc, id_var, '_FillValue', fill_val)
 pct_3D = ana_3Dpctvar  
 where (abs(ana_3Dpctvar) > 100)  pct_3D = fill_val   
 nc_status = nf90_put_var(lu_nc, id_var, pct_3D)
 if (nc_status .ne. 0) return
 
 !! copy a priori variance field in nc file 
 var_name = PARAM//'_ANO_ERR'
 nc_status = nf90_inq_varid(lu_nc, var_name, id_var)
 if (nc_status .ne. 0) return
 nc_status = nf90_get_att(lu_nc, id_var, 'scale_factor', scal_fact)
 nc_status = nf90_get_att(lu_nc, id_var, '_FillValue', fill_val)
 imax = abs(fill_val) - 1
 xmax = imax*scal_fact 
 std_3D = ana_3Dvarapr/scal_fact 
 where (abs(ana_3Dvarapr) > xmax)  std_3D = fill_val
 where (std_3D < 0 )  std_3D = fill_val
 nc_status = nf90_put_var(lu_nc, id_var, std_3D)
 if (nc_status .ne. 0) return
       
 
 nc_status = nf90_close (lu_nc)

return
end
