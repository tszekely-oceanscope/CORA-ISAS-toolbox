 
!===============================================================================
 
  subroutine  OA_ncwrite_dat (nc_fname, PARAM, nb_level, nb_data,  &
                     data_residl, data_errUR, nc_status )
 
!===============================================================================
!
!  Objet: Writes analysis results in nc file prepared by pre-OA
!  -----    
!
!  Version: 4.0 f95 - mpi
!  -------
! history:
! V4.0	: 22/05/2006: 	F. Gaillard	- F95-MPI
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
                 nb_data	! number of data
 
real*8  :: data_residl(nb_level,nb_data), & ! data table
           data_errUR(nb_level,nb_data)     ! error due to unresolved scales
				
!  OUTPUTS:
!  ------
integer*4 :: nc_status    ! status index after reading nc_file 


!  LOCALS:
!  ------
integer*4    :: lu_nc, io_mode, id_var

real*8       :: fill_val
		
character*20   :: var_name


!  DYNAMICALLY ALLOCATED ARRAYS

 real*4, allocatable :: data_tab(:,:)


!  --------------------------------
!  Program starts
!  --------------------------------

! Opens nc file 
 lu_nc     = 30    
 nc_status = nf90_open(nc_fname, nf90_write, lu_nc)
 write(*,*) nc_status
 if (nc_status .ne. 0) return
 
 allocate (data_tab(nb_level,nb_data))
 
 var_name = PARAM//'_RESI'
 nc_status = nf90_inq_varid(lu_nc, var_name, id_var)
 nc_status = nf90_get_att(lu_nc, id_var, '_FillValue', fill_val)
 data_tab = data_residl
 where(abs(data_tab) .GT. 900) data_tab = fill_val
 nc_status = nf90_put_var(lu_nc, id_var, data_tab)
 
 
 var_name = PARAM//'_ERUR'
 nc_status = nf90_inq_varid(lu_nc, var_name, id_var)
 nc_status = nf90_get_att(lu_nc, id_var, '_FillValue', fill_val)
 !write(*,*) id_var, fill_val
 data_tab = data_errUR
 where(abs(data_tab) .GT. 900) data_tab = fill_val
 nc_status = nf90_put_var(lu_nc, id_var, data_tab)

 
 nc_status = nf90_close (lu_nc)

 return
 end
