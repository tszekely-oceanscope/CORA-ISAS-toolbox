 
!===============================================================================
 
 subroutine OA_ncreaddata(nc_fname_dat, PARAM, jjul_est, qc_max, fill_val, &
                          nb_level, nb_data, nc_status, data_xyzt,         &
                          data_val, data_cli, data_cli_std, data_err)

!===============================================================================
!
!  Objet: Reads data in nc file prepared by pre-OA
!  -----    
!
!  Version: 4.0 f95 - mpi
!  -------
! history:
! V4.00	: 22/05/2006: 	F. Gaillard	- F95-
! V4.01	: 15/05/2007: 	F. Gaillard	- F95-
! V4.02	: 29/08/2007: 	F. Gaillard	- F95-
!         problems reading QC, this part is commented - QC are ignored
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
character*(*) :: nc_fname_dat 	! Name of nc file holding the data
character*4   :: PARAM 		! Type of variable (TEMP or PSAL)
real*8        :: jjul_est	! Estimate date (relative julian day)
real*8        :: fill_val	! 
integer*4     :: qc_max, &  	! max QC used
                 nb_level, & 	! Number of data levels 
                 nb_data   	! Number of data profiles 
         

!  OUTPUTS:
!  ------
integer*4 :: nc_status          ! status index after reading nc_file 
real*8   :: data_xyzt(nb_data,4)	! array of position and date of data
   				! xyt_dat(nb_level, 1): longitude
   				! xyt_dat(nb_level, 2): latitude	
				! xyt_dat(nb_level, 3): bottom depth (0 if not def)
   				! xyt_dat(nb_level, 4): julian day 
				! (relative to estim_date)
real*8   :: data_val(nb_level,nb_data)	! array of parameter value
real*8   :: data_cli(nb_level,nb_data)	! array of parameter reference climatology
real*8   :: data_cli_std(nb_level,nb_data)	! array of parameter reference climatology standard deviation
real*8   :: data_err(nb_level,nb_data)	! array of parameter measurement error


!  LOCALS:
!  ------
integer*4    :: lu_nc, io_mode, id_var, ln1
		
character*8   :: nom_var_lat
character*9   :: nom_var_lon
character*4   :: nom_var_jul
character*4   :: var_name
character*9   :: cli_name
character*13  :: cli_std_name
character*11  :: err_name 
character*7   :: qc_name


!  DYNAMICALLY ALLOCATED ARRAYS

real*4 ,   allocatable :: wk_sp4(:,:)
real*8,    allocatable :: wk_sp8(:)
integer*1, allocatable :: tab_qc(:,:)
integer*4, allocatable :: val_qc(:,:)


! tests dimensions
        character *16   :: name
        integer*4     :: xtype, ndims, ilen
        integer*4 :: dimids(2)
        integer*4 :: nAtts




!  --------------------------------
!  Program starts
!  --------------------------------


 nom_var_lat  = 'LATITUDE'
 nom_var_lon  = 'LONGITUDE'
 nom_var_jul  = 'JULD'
  
 var_name = PARAM 
 cli_name = PARAM//'_CLMN'
 cli_std_name = PARAM//'_CLSD'
 err_name = PARAM//'_ERME'
 qc_name  = PARAM//'_QC'

 ! Opens nc file 
lu_nc     = 30    
ln1 = index(nc_fname_dat, '.', back=.true.) + 2 
nc_status = nf90_open(nc_fname_dat(1:ln1), nf90_nowrite, lu_nc)
if (nc_status .ne. 0) return



!  Latitude, Longitude, date (jjul_rel):
!  -------------------------------------

 allocate(wk_sp8(nb_data))

 nc_status = nf90_inq_varid(lu_nc, nom_var_lon, id_var)
 nc_status = nf90_get_var(lu_nc, id_var, wk_sp8)
 if (nc_status .ne. 0) return
 data_xyzt(1:nb_data,1) = wk_sp8 

 nc_status = nf90_inq_varid(lu_nc, nom_var_lat, id_var)
 nc_status = nf90_get_var(lu_nc, id_var, wk_sp8)
 if (nc_status .ne. 0) return
 data_xyzt(1:nb_data,2) = wk_sp8 

 data_xyzt(1:nb_data,3) = 0

 nc_status = nf90_inq_varid(lu_nc, nom_var_jul, id_var)
 nc_status = nf90_get_var(lu_nc, id_var, wk_sp8)
 if (nc_status .ne. 0) return
 data_xyzt(1:nb_data,4) = wk_sp8 - jjul_est
 ! write (*,*) 'coord: ' , data_xyzt(1,:)

  deallocate(wk_sp8) 


! donnees et erreur
! -----------------
 nc_status = nf90_inq_varid(lu_nc, var_name, id_var)
 if (nc_status .ne. 0)  then
     nb_data = 0
     nc_status = 100
     return
 endif

 nc_status = nf90_inquire_variable(lu_nc, id_var, var_name, xtype, &
                                   ndims, dimids, nAtts)
 nc_status = nf90_inquire_dimension(lu_nc, dimids(1), name, len = ilen)
 if (ilen .ne. nb_level) then
    write (*,*) 'Dimension error: ', name, ilen, nb_level 
    return
 endif
 nc_status = nf90_inquire_dimension(lu_nc, dimids(2), name, len = ilen)
 if (ilen .ne. nb_data) then
    write (*,*) 'Dimension error: ', name, ilen, nb_data 
    return
 endif

 allocate(wk_sp4(nb_level, nb_data))
 allocate(tab_qc(nb_level, nb_data))
 allocate(val_qc(nb_level, nb_data))
 
 nc_status = nf90_inq_varid(lu_nc, var_name, id_var)
 nc_status = nf90_get_var(lu_nc, id_var, wk_sp4)

!write (*,*) var_name, id_var, nc_status

 if (nc_status .ne. 0)  then
    return
 endif
 data_val = dble(wk_sp4)

 nc_status = nf90_inq_varid(lu_nc, cli_name, id_var)
 nc_status = nf90_get_var(lu_nc, id_var, wk_sp4)
 data_cli = dble(wk_sp4) 

 nc_status = nf90_inq_varid(lu_nc, cli_std_name, id_var)
 nc_status = nf90_get_var(lu_nc, id_var, wk_sp4)
 data_cli_std = dble(wk_sp4) 

 nc_status = nf90_inq_varid(lu_nc, err_name, id_var)
 nc_status = nf90_get_var(lu_nc, id_var, wk_sp4)
 data_err = dble(wk_sp4) 
 
 nc_status = nf90_inq_varid(lu_nc, qc_name, id_var)
 nc_status = nf90_get_var(lu_nc, id_var, tab_qc)
 val_qc = tab_qc 

where (val_qc .GT. qc_max) data_err = fill_val
where (val_qc .GT. qc_max) data_val = fill_val
where (val_qc .GT. qc_max) data_cli = fill_val 

deallocate(wk_sp4, tab_qc, val_qc)    

nc_status = nf90_close (lu_nc)

!write (*,*) 'data_val (1,1:nb_data)'
!write (*,*) data_val(1,1:nb_data)
!write (*,*) 'data_val(1:nb_level,1)'
!write (*,*) data_val(1:nb_level,1)
!write (*,*) 'data_cli (1,1:nb_data)'
!write (*,*) data_cli(1,1:nb_data)
!write (*,*) 'data_cli(1:nb_level,1)'
!write (*,*) data_cli(1:nb_level,1)



return
end
