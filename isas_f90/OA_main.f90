
!  !-------------------------------------------------------------------------------
!
!
      PROGRAM OA_main
!
!  !-------------------------------------------------------------------------------
!
! history:
! V4.00	: 22/05/2006: 	F. Gaillard	- F95
!
! V4.01	: 14/05/2007: 	F. Gaillard	- F95
!
! V4.02	: 30/08/2007: 	F. Gaillard	- F95
!         cov_max added as external parameter
!
! V5.01	: 29/09/2009: 	F. Gaillard	- F95
!	allows filename <=160
!
! Calls: 
!	netcdf library
!	OA_anaarea
! !-------------------------------------------------------------------------------

implicit none

integer*4,  PARAMETER :: lurstd = 5
		
character*160 :: conf_fname, nc_dirname, log_dirname,           &
                 varapr_fname, covsca_fname, bathy_fname,       &
                 run_name, nc_fname, log_fname, err_fname, fname0
character*16 :: oa_version
character*4  :: PARAM
character*3  :: nam_area
integer*4    :: lu_conf, lu_err, ios, iarea, nb_area, qc_max, mx_std,    &  
                conf_xyzt(4), i, ii, ln1, ln2, ln3, ln4, ln5
real*8       :: var_weigh(3), cov_max, covar_ms_t, covar_ls(3), &
                tps_tot, fact, tps_all


 oa_version = 'OA-Version 7.0'
 tps_all = 0
 
 read(lurstd,'(a)') nc_dirname 
 read(lurstd,'(a)') log_dirname  
 read(lurstd,'(a)') conf_fname 
 ln1 = index(conf_fname, '$') - 1 
 write (*,*) 
 write (*,*) 'conf_fname: ', conf_fname(1:ln1)


! Open config file
! ----------------
 lu_conf = 22
 open(unit=lu_conf, file= conf_fname(1:ln1), status='old', form='formatted',  & 
                                  access='sequential',iostat=ios)
 if (ios  .ne. 0) then
    write(*,*) ' OA_main: error opening config file, ios = ', ios 
    stop
 endif
 
   read(lu_conf,'(a4)') PARAM 
   read(lu_conf,'(a)') varapr_fname   
   read(lu_conf,'(a)') covsca_fname      
   read(lu_conf,'(a)') bathy_fname   
   read(lu_conf, *) (covar_ls(i),  i = 1, 3)
   read(lu_conf, *) covar_ms_t 
   read(lu_conf, *) (var_weigh(i), i = 1, 3)    
   read(lu_conf, *) (conf_xyzt(i), i = 1, 4)    
   read(lu_conf, *) fact 
   read(lu_conf, *) qc_max, mx_std
   read(lu_conf, *) cov_max
 close(lu_conf)

 covar_ls(1:2) = covar_ls(1:2)*1000	! Converts to meters

 read(lurstd, *) nb_area  

 do ii = 1, nb_area

    read(lurstd,'(a)')  fname0 
    write (*,*) fname0
    ln1 = index(fname0, '.', back=.true.) - 1   
    ln2 = index(fname0, ' ') - 1   
    ln1 = min(ln1,ln2)
    nam_area = fname0(ln1-11:ln1-9)
    read(nam_area, '(i3)') iarea
    
    ln3 = index(nc_dirname, '$') - 1 
    ln4 = index(log_dirname, '$') - 1  

    if (ii .eq. 1) then

!  Opens error file 
!  ---------------  
        lu_err = 21    
        run_name = fname0(1:ln1-12)//fname0(ln1-3:ln1)
        ln5 = len_trim(run_name)
 
        err_fname = 'err/'//run_name(1:ln5)//'.err'
        open(unit=lu_err, file=err_fname, status='unknown', form='formatted',  & 
                                   access='sequential',iostat=ios)
        if (ios  .ne. 0) then
            write(*,*) ' OA_main: Error opening error file, ios = ', ios 
            stop	
        endif
        
        write(*,10110) run_name(1:ln5), oa_version, PARAM, nb_area

        write(lu_err,10110) run_name(1:ln5), oa_version, PARAM, nb_area
        write(lu_err,10120) varapr_fname(1:len_trim(varapr_fname)), &
                            covsca_fname(1:len_trim(covsca_fname)), &
                            bathy_fname(1:len_trim(bathy_fname))
        write(lu_err,10130) (covar_ls(i),  i = 1, 3), covar_ms_t, &
                            (var_weigh(i), i = 1, 3), &
                            (conf_xyzt(i), i = 1, 4), &   
                            fact, qc_max, mx_std,  &
                            cov_max


  endif


  ! Computes analyzed field
  ! -----------------------
    log_fname = log_dirname(1:ln4)//fname0(1:ln1)//'.log'  !  
    nc_fname  = nc_dirname(1:ln3)//fname0 (1:ln1)//'.nc'
    write(lu_err,10140) fname0(1:ln1)
    call OA_anaarea (lu_err, log_fname, nc_fname,                        &
                     varapr_fname, covsca_fname, bathy_fname,            &
		     PARAM, iarea, mx_std, fact, qc_max, var_weigh,      &
		     cov_max, conf_xyzt, covar_ls, covar_ms_t, tps_tot)

    tps_all = tps_all + tps_tot
 
 enddo

 
write (*,*) tps_all 
write(lu_err,*) tps_all
close(lu_err)

 stop	


10110  format (/ , ' ==============================================================' /  &
                   '     Run: ', a ,   / '           ' a                             /, & 
                   ' ==============================================================' // &
                   '    PARAM: ', a, ', nb_area: ', i4 /)

10120  format (/, '  Configuration files used: ' / &               
                  '  ------------------------' / &
               '  A-priori variance :' /, a /&
               '  Covariance scales :' /, a /&
               '  Bathymetry :' /, a /)

10130  format (/, ' Parameters for optimal estimation: '/ &
               ' ---------------------------------' / &
               ' Large scale covariance (x,y,t): ', 2f10.1, f8.1 / &
               ' Meso scale covariance (t)     : ', f8.1 / &
               ' Variance weights (LS, MS, UR) : ', 3f6.1 / &
               ' Covariances used (x, y, z, t) : ', 4i2 / &
               ' Factor multiplying apr-var    : ', f4.1 / &
               ' qc_max, mx_std:               ', 2i4, / &
               ' Oversampling: cov_max', f8.5, /)

10140  	 format (/, '  Analysis name:', a)

 end PROGRAM OA_main
