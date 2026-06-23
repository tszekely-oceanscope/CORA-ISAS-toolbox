
!
!==============================================================================
! 
 subroutine OA_anaarea (lu_err, log_fname, nc_fname,                        &
                        varapr_fname, covsca_fname, bathy_fname,            &
			PARAM, iarea, mx_std, fact, qc_max, var_weigh,      &
			cov_max, conf_xyzt, covar_ls, covar_ms_t, &
                        tps_tot)
!
!==============================================================================
!
!  Objet: Analysis of a field on a regular grid
!  -----
!
! history:
!--------
! V1	: 1999 - F. Gaillard, C. Lagadec - creation matlab
! V2	: 2002 - E.Autret, F. Gaillard - matlab
! V3	: 2004 - E.Autret - matlab
! V4.00	: 07/06/2006: 	F. Gaillard	- F95
! V4.01	: 09/05/2007: 	F. Gaillard	
!         uses variances and scales defined at each analysis point
! V4.02	: 31/08/2007: 	F. Gaillard	
!	  corrected bug on data covariances
!
! Calls: 	
!       OA_ncreaddim
!       OA_ncreadstat
!	OA_ncreaddata
!	OA_covini
!	OA_calsol
!	OA_ncwrite_fld, 
!       OA_ncwrite_dat
!
!  Reads and writes in [fname'.nc']
!  Writes in [fname '.log']		
!
!  Method:
!  =======
!  X_ana = X_0 + Cmd (Cdd + Co)-1 (Y _ X_0)
! Cdd(i,j) = Sig_i*Sig_j*F(X_i,X_j)
! Sig_i*Sig_i = var_ana_apr_i
!
!   Computing apriori variances 
!   ---------------------------
!  Var_tot = Var_LS + Var_MS + Var_UR + Var_ME
!   	Var_LS = large scale
!	Var_MS = mesoscale
! 	Var_UR = subgrid or unresolved scales
! 	Var_ME = measurement error
!
! Var_tot is obtained from statistics on data base !
! Var_stat = var_tot - Var_ME, is stored in area nc file 
! Var_stat = Wgh1*Var_stat + Wgh2*Var_stat + Wgh3*Var_stat 
! sum(Wgh) = 1, Wgh are defined in conf_xyzt file
!
! var_ana_apr  = Var_LS + Var_MS
! var_data_apr = Var_UR + Var_ME
!
! Var_LS = Wgh1*Var_stat	needed on analysis grid
! Var_MS = Wgh2*Var_stat    	needed on analysis grid
!
! Var_UR = mean(Var_0*Wgh3) 	needed at data points
! Var_ME : read in data file	needed at data points
!
!==============================================================================


implicit none

integer*4,  PARAMETER :: mx_gauss = 2


!  INPUTS:
!  ------
character*(*) :: nc_fname           ! Name of area file (used for nc (R/W) 
character*(*) :: log_fname          ! Name of log file (W)
character*(*) :: varapr_fname, covsca_fname, bathy_fname
character*4   :: PARAM 	            ! Type of variable (TEMP or PSAL)
integer*4    :: lu_err, &
                iarea , &	    ! index of area to be processed 
                qc_max 
real*8       :: var_weigh(3), & ! (relative weighr of LS, MS, SG)
                cov_max,      & ! covariance max (calcul d'oversampling)
                covar_ls(3), covar_ms_t, &
		fact 
integer*4    :: conf_xyzt(4) 	! defines covariance computation 
				! (switches for x, y, z, t)
			     	! 1 1 1 1: x y z t on

!  OUTPUTS:
!  -------
real*8       :: tps_tot		! CPU time for area processing


!  LOCALS:
! -------
integer*4    :: nb_dat_min,        &  ! Min number of data needed for analysis
                n_gauss,   mx_std, & ! Number of gaussian, criteria for std max
		lu_log, lu_nc,     &  ! Logical unit for files
		ios, i_data, i_level, nc_status, ifile, &
		nb_data, nb_level, mx_ana, nx_ana, ny_ana, nb_ana, &
		nx_bathy, ny_bathy,    &
		nbok_data, nbok_data2, nbok_ana, nb_elim, &
		i, j, ilev, ii, jj, ij, id, im, ijok, i_err, &
		i_test, ln1, ln2
		
real*8  :: fill_val, var_min, &    
           jjul_est,  wght_rel(3), w12, &	 
	   var_ij(2), ano_max, cond_i, dat_min, dat_max, &
           err_max, err_min, err_max_tot, err_min_tot
	   
real*4    tps_zero, tps_beg, tps_end, tps_ini, tps_cal 

character*(80) :: nc_fname_fld, nc_fname_dat ! Name of nc file holding the data

!logical :: 


!  DYNAMICALLY ALLOCATED ARRAYS
!  ----------------------------
! Tables defined at analysis gridpoints (2D)
! apriori info:
real*8 , allocatable :: lat_ana(:), lon_ana(:), covar_ms(:,:,:), bath_ana(:,:), &
                        ana_3Dvartot(:,:,:), dxx(:), dyy(:)
real*8 , allocatable :: lat_bat(:), lon_bat(:), covar_all(:,:,:), bath_all(:,:)
integer*4,  allocatable :: msk_area(:,:), i0(:), j0(:)

! results:
real*8 , allocatable :: ana_3Dfld(:,:,:), ana_3Dpctvar(:,:,:), &
                        ana_3Dvarapr(:,:,:)


! tables defined at data points
real*8 , allocatable 	:: data_xyzt(:,:),   data_val(:,:),  &
                           data_cli(:,:), data_cli_std(:,:), &
			   data_err(:,:), data_fldvar(:,:),  &
			   data_covms(:,:), data_errUR(:,:),   &
			   data_resid(:,:)

! Covariance matrices
real*8 , allocatable   	:: ana_xyz(:,:), ana_fldvar(:,:), ana_covms(:,:)
real*8 , allocatable   	:: F_dd(:,:,:), F_md(:,:,:)
integer*4, allocatable 	:: idx_ana(:,:), isok(:,:)


real*8 , allocatable 	:: dx(:), dy(:), dd(:)

real*8 , allocatable   	::  std_lev(:), an_vpri(:), an_fldi(:), &
                            an_sigpri(:,:), an_vpsi(:)
 			
real*8 , allocatable   	:: data_ano_i(:), data_var_i(:), dat_sigfldi(:,:), &
                           var_dat_UR_i(:), data_res_i(:)

real*8 , allocatable   	::  C_dd_i(:,:),  C_md_i(:,:)
	
integer*4, allocatable 	:: iok_data(:), iok_data2(:), iok_ana(:)



!  --------------------------------
!  Program starts
!  --------------------------------
 
 i_test = 0 ! i_test = 1 for more outputs in log file


 nb_dat_min = 3
 n_gauss = mx_gauss
 fill_val = 99999.0D0
 
 var_min = 0.0001



 call cpu_time(tps_zero)


! Open log file 
!  -------------  
 lu_log = 20    
 open(unit=lu_log, file=log_fname, status='unknown', form='formatted',  & 
                                   access='sequential',iostat=ios)
 if (ios  .ne. 0) then
    write(lu_err,10902) log_fname, ios 
    write(*,10902) log_fname, ios 
   return	
 endif



! Reads dimensions
! ---------------- 
 ln1 = index(nc_fname, '.') - 1
 ln2 =  ln1 - 8   
 nc_fname_dat = nc_fname(1:ln2)//'dat_'//nc_fname(ln2+5:ln1)//'.nc'
 nc_fname_fld = nc_fname(1:ln2)//'fld_'//nc_fname(ln2+5:ln1)//'.nc'
 nc_status = 0
 call OA_ncreaddim(nc_fname_fld, nc_fname_dat, bathy_fname, &  
                   nc_status, nb_data, nx_ana, ny_ana, nb_level, &
                   jjul_est, nx_bathy, ny_bathy)
 if (nc_status  .ne. 0) then      
    ifile = int(nc_status/100)
    nc_status = nc_status - 100*ifile
    write(*,10910) ifile, nc_status
    write(lu_err,10910) ifile, nc_status
    write(lu_log,10910) ifile, nc_status
    close (lu_log)
    return
 endif
 

 write (*, 10100) iarea, nb_data, nb_level, nx_ana, ny_ana, &
                  nx_bathy, ny_bathy 
 write (lu_log, 10100) iarea, nb_data, nb_level, nx_ana, ny_ana, &
                       nx_bathy, ny_bathy

 if (nb_data < nb_dat_min) then 
     write (*, 10105) nb_data, nb_dat_min
     write (lu_log, 10105) nb_data, nb_dat_min
     close(lu_log)
     return
 end if

 
 ! Reads variances and covariances scales 
 ! ---------------------------------------
 
 allocate (lat_ana(ny_ana), lon_ana(nx_ana))
 allocate (covar_ms(nx_ana,ny_ana,2), msk_area(nx_ana,ny_ana), &
           bath_ana(nx_ana,ny_ana))

 allocate(ana_3Dvartot(nx_ana, ny_ana, nb_level), &   
          ana_3Dfld(nx_ana, ny_ana, nb_level),    &   
          ana_3Dpctvar(nx_ana, ny_ana, nb_level), &   
          ana_3Dvarapr(nx_ana, ny_ana, nb_level))

 allocate (lat_bat(ny_bathy), lon_bat(nx_bathy),   &
           covar_all(nx_bathy,ny_bathy,2), bath_all(nx_bathy,ny_bathy))


 nc_status = 0
 call OA_ncreadstat(nc_fname_fld, varapr_fname, covsca_fname, bathy_fname, &
			     PARAM, fact, iarea, nx_ana, ny_ana, nb_level, &
                             nx_bathy, ny_bathy, &
                             nc_status, lat_ana, lon_ana, ana_3Dvartot,    &
			     covar_ms, bath_ana, msk_area, &
                             lat_bat, lon_bat, covar_all, bath_all)


 if (nc_status  .ne. 0) then      
    ifile = int(nc_status/100)
    nc_status = nc_status - 100*ifile
    write(*,10910) ifile, nc_status
    write(lu_err,10910) ifile, nc_status
    write(lu_log,10910) ifile, nc_status
    close (lu_log)
    return

 else
    write (*, *) '   OA_anaarea: statistics OK '   
    write (lu_log, *) '   OA_anaarea: statistics OK '
 endif

 ana_3Dvarapr(1:nx_ana, 1:ny_ana, 1:nb_level) = fill_val
 ana_3Dfld(1:nx_ana, 1:ny_ana, 1:nb_level)    = fill_val
 ana_3Dpctvar(1:nx_ana, 1:ny_ana, 1:nb_level) = fill_val
 covar_ms(1:nx_ana,1:ny_ana,1:2) = covar_ms(1:nx_ana,1:ny_ana,1:2)*1000	! (km--> meter)

 
 if (i_test .EQ. 1) then
     i = 5
     write(lu_log,'(/)')
     write(lu_log,*) 'i=5 msk_area bath_ana  cov_ms_X cov_ms_Y ana_3Dvartot (j = ny_ana)'
     do j = 1, ny_ana
        write (lu_log,10801) j, msk_area(i,j), bath_ana(i,j),& 
                               covar_ms(i,j,1), covar_ms(i,j,2), ana_3Dvartot(i,j,1)
     enddo
 endif


 ! Read data
 ! ---------
 allocate(data_xyzt(nb_data,4),            &
          data_val(nb_level,nb_data),      &
          data_cli(nb_level,nb_data),      &
          data_cli_std(nb_level,nb_data),  &
	  data_err(nb_level,nb_data),      &
	  data_resid(nb_level, nb_data),   &  
          data_errUR(nb_level, nb_data),   &
	  data_fldvar(nb_level,nb_data),   &
 	  data_covms(nb_data,2))
 allocate(dyy(ny_bathy), dxx(nx_bathy), j0(ny_bathy), i0(nx_bathy))
	
 call OA_ncreaddata(nc_fname_dat, PARAM, jjul_est, qc_max, fill_val, &
                          nb_level, nb_data, nc_status, data_xyzt,   &
                          data_val, data_cli, data_cli_std, data_err)

	    
 if (nc_status  .ne. 0) then    
    write (*, 10912) nc_status     
    write (lu_err, 10912) nc_status 
    write (lu_log, 10912) nc_status 
    close (lu_log)
    return
    
 else
    write (*, *) '   OA_anaarea: data OK '   
    write (lu_log, *) '   OA_anaarea: data OK '
 endif

 
 do ii = 1, nb_data
     dxx = abs(lon_bat - data_xyzt(ii,1))
     i0 = minloc(dxx)   ! find nearest grid point (i,j)
     dyy = abs(lat_bat - data_xyzt(ii,2))
     j0 = minloc(dyy) 
     data_fldvar(1:nb_level,ii) =  &
                  fact*data_cli_std(1:nb_level,ii)*data_cli_std(1:nb_level,ii)
     data_xyzt(ii,3) = bath_all(i0(1),j0(1))
     data_covms(ii,1:2)  = covar_all(i0(1),j0(1),1:2)*1000	! (km--> meter)
 enddo

 
 data_resid(1:nb_level, 1:nb_data)  = fill_val
 data_errUR(1:nb_level, 1:nb_data)  = fill_val


 deallocate (dyy, dxx, i0, j0)

 deallocate (lon_bat, lat_bat, bath_all, covar_all)


 if (i_test .EQ. 1) then
      write(lu_log,'(/)')
      write (lu_log,*) ' Checking data:'
     do i = 1, nb_data     
        write (lu_log,10802) i, data_xyzt(i,1:3), data_covms(i,1:2)
        do j = 1,nb_level,20
            write (lu_log,10803) j, data_val(j,i), data_err(j,i), &
                                            data_cli(j,i), data_cli_std(j,i)         
        enddo	
     enddo	
 endif  


 !  sort analysis points in 1D 
 !  -------------------------
 mx_ana = nx_ana*ny_ana
 allocate(ana_xyz(mx_ana,3), idx_ana(mx_ana,2), &
          ana_fldvar(nb_level,mx_ana), ana_covms(mx_ana,2))
 ii = 0
 do j = 1, ny_ana
    do i = 1, nx_ana
       if (msk_area(i,j) .EQ. iarea) then
          ii = ii + 1
          ana_xyz(ii,1) = lon_ana(i)         
	  ana_xyz(ii,2) = lat_ana(j)         
	  ana_xyz(ii,3) = bath_ana(i,j) 
	  ana_fldvar(1:nb_level,ii) = ana_3Dvartot(i,j,1:nb_level)
	  ana_covms(ii,1) = covar_ms(i,j,1)
	  ana_covms(ii,2) = covar_ms(i,j,2)
          idx_ana(ii,1) = i
	  idx_ana(ii,2) = j        
       endif
    enddo
 enddo
 nb_ana = ii

 deallocate (covar_ms, ana_3Dvartot)


 if (i_test .EQ. 1) then	
     write(lu_log,*) ' Analysis: mx_ana, nb_ana : ', mx_ana, nb_ana 
    
     write (lu_log,*) 'lon_ana' 
     write (lu_log,*)  lon_ana 
     write (lu_log,*) 'lat_ana'
     write (lu_log,*)  lat_ana

     write(lu_log,*) ' ana_xyz(i,:), ana_covms(i,1:2)'
     do i = 1, nb_ana, 100
        write(lu_log,10802) i, ana_xyz(i,1:3), ana_covms(i,1:2)
     enddo 
     write(lu_log,*)  'ana_fldvar(1,1:nb_ana)'
     write(lu_log,'(12f6.2)') ana_fldvar(1,1:nb_ana)
 endif
  


! Prepares covariance matrices (data-data and ana-data):
! ------------------------------------------------------
 
 allocate (F_dd(nb_data,nb_data,mx_gauss), F_md(nb_ana,nb_data,mx_gauss))

 i_data = nb_data
 
 call cpu_time(tps_beg)	   
 call OA_covini (data_xyzt, i_data, nb_data, data_covms, &
                 ana_xyz, mx_ana, nb_ana, ana_covms, &
                 n_gauss, conf_xyzt, covar_ls, covar_ms_t, &
		 F_dd, F_md)
 call cpu_time(tps_end)
 tps_ini = tps_end - tps_beg
 write (lu_log, 10110) tps_ini

 deallocate (ana_covms, data_covms)  


 if (i_test .EQ. 1) then
      i = 1  
      write(lu_log,*) nb_ana, i       
      write(lu_log,*)  'F_dd(1:nb_data,1,1)'
      write(lu_log,10804) F_dd(1:nb_data,1,1)               
      write(lu_log,*) 'F_md(i,1:nb_data,1)'  
      write(lu_log,10804) F_md(i,1:nb_data,1)               
      write(lu_log,*)  'F_dd(1:nb_data,1,2)'
      write(lu_log,10804) F_dd(1:nb_data,1,2)               
      write(lu_log,*) 'F_md(i,1:nb_data,2)'  
      write(lu_log,10804) F_md(i,1:nb_data,2)                             
 endif


!
! ============================================================================
!                           Loop over analysis levels
! ============================================================================ 

 allocate (iok_data(nb_data), iok_data2(nb_data), iok_ana(nb_ana))
 allocate (std_lev(nb_ana)) 
 allocate (an_vpri(nb_ana), an_fldi(nb_ana), &
           an_sigpri(nb_ana,2), an_vpsi(nb_ana)) 
 allocate (data_ano_i(nb_data), data_var_i(nb_data), &
           dat_sigfldi(nb_data,2), var_dat_UR_i(nb_data), &
           data_res_i(nb_data))     
 allocate (C_dd_i(nb_data,nb_data), C_md_i(nb_ana,nb_data))    
 
 wght_rel = var_weigh/sum(var_weigh)
 w12 = var_weigh(1) + var_weigh(2)


 do ilev = 1,nb_level 
 !do ilev = 2,2 
 ! looks for valid analysis points     
    i = 0 
    do j = 1, nb_ana    
      if (ana_fldvar(ilev,j) .GT. var_min) then
              i = i + 1
	      iok_ana(i) = j
       endif  
    enddo
    nbok_ana = i
 
    if (nbok_ana .GT. 0) then
       std_lev(1:nbok_ana) = sqrt(ana_fldvar(ilev,iok_ana(1:nbok_ana)))  
       ano_max = mx_std*maxval(std_lev(1:nbok_ana))
 
! looks for valid data 
        i = 0!!! remplace un find  (on peut surement faire mieux) 
        do j = 1, nb_data   
           if (abs(data_val(ilev,j)) .LT. 9000 .AND. &
               data_fldvar(ilev,j) .GT. var_min) then
!write(*,*) data_val(ilev,j), data_cli(ilev,j)
               if (abs(data_val(ilev,j) - data_cli(ilev,j)) .LT.  ano_max) then
                  i = i + 1
	          iok_data(i) = j
               endif
           endif  
        enddo
        nbok_data = i
    endif

    write (lu_log, 10120) ilev, nbok_ana, nbok_data


!  Identifies oversampling data and removes them
    call OA_oversamp (F_dd, nb_data,mx_gauss, var_weigh,w12, cov_max, &
                          iok_data, nbok_data, iok_data2, nbok_data2) 
    write (lu_log, 10121) nbok_data2
    nbok_data = nbok_data2
    iok_data = iok_data2
  

    if (nbok_data < nb_dat_min) then  ! Checks data availability
            write (*, 10122) ilev, nbok_data, nb_dat_min       
            write (lu_log, 10122) ilev, nbok_data, nb_dat_min  
    else  
      
        an_sigpri(1:nbok_ana,1) = sqrt(wght_rel(1))*std_lev(1:nbok_ana)   
        an_sigpri(1:nbok_ana,2) = sqrt(wght_rel(2))*std_lev(1:nbok_ana)
        an_vpri(1:nbok_ana) = an_sigpri(1:nbok_ana,1)*an_sigpri(1:nbok_ana,1) &
                            + an_sigpri(1:nbok_ana,2)*an_sigpri(1:nbok_ana,2)
			 
        do i = 1, nbok_data
           ii = iok_data(i)
		 
! extract valid data
           data_ano_i(i) = data_val(ilev,ii) - data_cli(ilev,ii)

 ! a priori variance on data    
           var_dat_UR_i(i) = wght_rel(3)*data_fldvar(ilev,ii)  
           data_var_i(i) = data_err(ilev,ii)*data_err(ilev,ii) + var_dat_UR_i(i)
       
! a priori variance on field at data points
	   dat_sigfldi(i,1:2) = sqrt(wght_rel(1:2)*data_fldvar(ilev,ii))
  
	   do j = 1, i	! triangle superieur
	      jj = iok_data(j)
	      var_ij = dat_sigfldi(j,:)*dat_sigfldi(i,:)
	      C_dd_i(j,i) = var_ij(1)*F_dd(jj,ii,1) + var_ij(2)*F_dd(jj,ii,2)
	   enddo

	   do j = 1, nbok_ana
	      jj = iok_ana(j)
	      var_ij = dat_sigfldi(i,:)*an_sigpri(j,:)
	      C_md_i(j,i) = var_ij(1)*F_md(jj,ii,1) + var_ij(2)*F_md(jj,ii,2)
	   enddo
!	write (*,*)  i, data_ano_i(i), data_var_i(i), C_dd_i(i,i)
 
        enddo   !	do i = 1, nbok_data
    
        dat_max = maxval(data_ano_i(1:nbok_data))
        dat_min = minval(data_ano_i(1:nbok_data))
       
        write (lu_log, 10125) ano_max, dat_min, dat_max

 
   ! Symetrise la matrice de covariance  C_dd   (pas nécessaire)
        do j = 1, nbok_data	 
           do i = j+1, nbok_data
              C_dd_i(i,j) = C_dd_i(j,i);	 
           enddo
        enddo

  
!  Computes analyzed field
!  -----------------------
        call cpu_time(tps_beg)

        id = nb_data
        im = nb_ana
        call OA_calsol(C_md_i, im, C_dd_i, id, data_var_i, an_vpri, &
                        data_ano_i, nbok_data, nbok_ana, &
		        an_fldi, an_vpsi, data_res_i, cond_i, i_err)   
        call cpu_time(tps_end)
        tps_cal = tps_end - tps_beg

        if (i_test .EQ. 1) then	    
            write(lu_log,*)  'an_fldi'
            write(lu_log,'(12f6.2)') an_fldi               
            write(lu_log,*)  'an_vpri'
            write(lu_log,'(12f6.2)') an_vpri   
            write(lu_log,*)  'an_vpsi'
            write(lu_log,'(12f6.2)') an_vpsi
            write(lu_log,*)  'an_vpsi/an_vpri'
            write(lu_log,'(12f6.2)') an_vpsi/an_vpri
        endif

      
        if (i_err .NE. 0) then
            write(*,10128)  ilev
            write (lu_err, 10128) ilev 
            write (lu_log, 10128) ilev
     
        else
           dat_max = maxval(an_fldi(1:nbok_ana))
           dat_min = minval(an_fldi(1:nbok_ana))

          do i = 1, nbok_data
             ii = iok_data(i)	    
             data_resid(ilev,ii)  = data_res_i(i)
             data_errUR (ilev,ii) = sqrt(var_dat_UR_i(i))  
          enddo 
	    
          do ij = 1,nbok_ana
             ii = iok_ana(ij)    
             i = idx_ana(ii,1)
             j = idx_ana(ii,2)
             ana_3Dfld(i,j,ilev)    = an_fldi(ij)       
	     ana_3Dvarapr(i,j,ilev) = an_vpri(ij)
             ana_3Dpctvar(i,j,ilev) =  100*an_vpsi(ij)/an_vpri(ij)
          enddo


	  err_max = maxval(ana_3Dpctvar(1:nx_ana,1:ny_ana,ilev),  &
                  mask = ana_3Dpctvar(1:nx_ana,1:ny_ana,ilev)< 200)
	  err_min = minval(ana_3Dpctvar(1:nx_ana,1:ny_ana,ilev))

          write (lu_log, 10126) dat_min, dat_max, err_min, err_max
          write (lu_log, 10130) cond_i, tps_cal
          if (err_min .le. 0) then
             write (lu_err, 10129) ilev, err_min, cond_i
          endif
 	


       endif	!	if i_err .NE. 0 then

    endif ! if (nbok_data < nb_dat_min)  
 enddo   ! do ilev = 1,nb_level
 err_max_tot = maxval(ana_3Dpctvar(1:nx_ana,1:ny_ana,1:nb_level), &
               mask = ana_3Dpctvar(1:nx_ana,1:ny_ana,1:nb_level)< 200)
 err_min_tot = minval(ana_3Dpctvar(1:nx_ana,1:ny_ana,1:nb_level)) 
 deallocate (iok_data, iok_ana)
 deallocate (std_lev) 
 deallocate (an_vpri, an_fldi, an_sigpri, an_vpsi) 
 deallocate (data_ano_i, data_var_i, dat_sigfldi, var_dat_UR_i, data_res_i)     
 deallocate (C_dd_i, C_md_i)
 
 if (i_test .EQ. 1) then
     write(lu_log,*) 'ana_3Dfld ana_3Dpctvar(1,j,1)'
     do j = 1, ny_ana
        write (lu_log,*) j, ana_3Dfld(5,j,1), ana_3Dpctvar(5,j,1)
     enddo
 endif 
      

! ============================================================================
!                       End of loop over analysis levels
! ============================================================================ 
    
!    Saves results in nc_file
!    ------------------------
 call OA_ncwrite_fld (nc_fname_fld, PARAM, nb_level, ny_ana, nx_ana, &
                      ana_3Dfld, ana_3Dpctvar, ana_3Dvarapr, nc_status )
 
 call OA_ncwrite_dat (nc_fname_dat, PARAM, nb_level, nb_data, &
                     data_resid, data_errUR, nc_status )
  
 deallocate (ana_3Dfld, ana_3Dpctvar, ana_3Dvarapr)
 deallocate (data_resid, data_errUR)  

 call cpu_time(tps_end)
 tps_tot = tps_end - tps_zero
 
 write (lu_log, 10140) err_min_tot, err_max_tot, tps_tot

 if (err_min_tot<0)  then 
   write (*, 10140) err_min_tot, err_max_tot, tps_tot  
   write (lu_err, 10140) err_min_tot, err_max_tot, tps_tot 
 endif

 close (lu_log)

 return


!  ------------
!   Messages
!  ------------

10100 	format (//  '*****   Area: ' i4 ',       Nb_profiles: ' i5 /  &
                    '                             Nb_level: ' i5 /    &
		    '      Nb_analysis points: (nlon,nlat): ' i5, i5/ &
                    '      Nb_Bathy points:    (nlon,nlat): ' i5, i5/) 
10105   format (    ' area not processed, Nb_data: ' i5  ' < ' i5 )
10110 	format (    '  cpu distance calculations:' f7.3)
10120   format (/   ' Level: ' i4 ', Nb_ana_points: ' i5 ',   Nb_data: ' i5)
10121   format (    '                                   Nb_ok_data: ' i5 )
10122   format (    ' Level: ' i4 ' not processed, Nb_data: ' i5  ' < ' i5 )
10125   format (    '     ano_max: ' f7.3   &
                    ', inov min: ' f7.3 ', inov max: ' f7.3 )
10126   format (    '                     ' &
                    '  fld  min: ' f7.3 ', fld  max: ' f7.3 / &
                    '                     ' &
                    '   err_min: ' f7.3 ',  err max: ' f7.3 )
10128   format (    '  Level: ' i4 ', Problem in Cal-sol')
10129   format (    '  Level: ' i4 ',  err_min: ' f7.3 ',   cond # ' e12.4 )
10130   format (    '     cond # ' e12.4 ',     cpu Analysis:' f7.3)
10140   format ( /   '   *****   End area: ', / &
                     '     err_min_tot: ' f7.3 ', err max_tot: ' f7.3 / &
                     '  cpu total area :' f9.3 /&
                     '  ********************'//) 
10801  format (2i6, f6.0, 3f12.1) 
10802  format (i5, 2f10.3, 3f10.0)
10803  format (i5, 4f8.3)
10804  format (12e10.2)

10902  format (' OA_anaarea: Error opening log file: ' / a, /'  Error code: ' i3)
10910  format (' Error in input files' / &
               '  0: data, 1: field, 2: bathymetry, 3: cov-scales, 4: variance'/ &
                 '   pb file file: ', i2 ', nc_status = ', i4 /)
10912  format ( ' Error in OA_anaarea, data nc file nc_status = ' i4)

 
end
