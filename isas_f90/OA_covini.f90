
!===============================================================================
  
  subroutine OA_covini (data_xyzt, i_data, nb_data, data_covms, &
                        ana_xyz, i_ana, nb_ana, ana_covms, &
                        n_gauss, config, covar_ls, covar_ms_t, &
		        F_dd, F_md)

!===============================================================================
!
!  Objet: Prepares covariance matrices by computing the distance dependent part
!  -----    
!
!
! V1	: 1999 - F. Gaillard, C. Lagadec - creation matlab
! V2	: 2002 - E.Autret, F. Gaillard - matlab
! V3	: 2003 - E.Autret - matlab
! V4.0	: 22/05/2006: 	F. Gaillard	- F95-
! V4.01	: 14/05/2007: 	F. Gaillard	- F95-
! V5.02	: 15/10/2009: 	F. Gaillard	- F95-
!        Correction for 180 longitude line
!!
!
!  Method:
!  =======
!  F(x_i,X_j) = sum_ig = 1, n_gauss (exp[- 1/2(D_ig**2]).
!  D_ig = normalized distance
!  D_ig**2 = dx**2/Lx_ig**2  + dy**2/Ly_ig**2  + dt**2/Lt_ig**2  + dh/Lh_ig**2 
!

!  INPUT PARAMETERS
!       data_xyzt(i_data ,4):positions and date of data points 
!		- data_xyzt(nb_data,1) = x_long en degrees, 
!		- data_xyzt(nb_data,2) = y_lati en degrees, 
!		- data_xyzt(nb_data,3) = bottom depth
!		- data_xyzt(nb_data,4) = date in day relatives to estim day
!       i_data : leading dimension of data_xyzt
!       nb_data : number of data
!       data_covms : mesoscale covariances at data points
!                      
!      ana_xyz(i_ana,3):positions et depth of analyzed points
!       i_ana : leading dimension of ana_xyzt
!       nb_ana : number of analyzed points
!       ana_covms: mesoscale covariances at analysis points
!
!       n_gauss : number of gaussian to sum up
!       conf_xyzt: swithes for covariance dependency (1 = on, 0 = off
!                  x, y, z, t)
!
!		covar_ls(1,i_gauss) = Lx_ig: covariance scale in x (km)
!		covar_ls(2,i_gauss) = Ly_ig: covariance scale in y (km)
!		covar_ls(3,i_gauss) = Ly_ig: covariance scale in t  (days)
!           covar_ms_t: time scale for mesoscale covariance
!
!
!  OUPUT PARAMETERS         
!       F_dd(i_data,i_data,n_gauss):
!		    covariance data-data for each i_gauss
!	F_md(i_ana,i_data,n_gauss):
!		    covariance model-data for each i_gauss
!
!  NOTES:
!
!  The amplitudes sigma_i^2 are taken into account in  AO_CAL.
!
!===============================================================================
 
 implicit none
 
!  INPUTS:
!  ------
integer*4    :: i_data, nb_data, i_ana, nb_ana, i_cov, n_gauss, config(4)
real*8       :: data_xyzt(i_data ,4), ana_xyz(i_ana,3), &
                data_covms(nb_data,2), ana_covms(i_ana,2), covar_ls( 3), covar_ms_t

!  OUTPUTS:
!  -------
real*8       :: F_dd(i_data,i_data,n_gauss), F_md(nb_ana,i_data,n_gauss)	


!  LOCALS:
! -------
integer*4    :: i, j, i_gauss 
		
real*8  :: pi19, Rterre, radeg, degrad, convR, l_H, xbid
	   

!logical :: 


!  DYNAMICALLY ALLOCATED ARRAYS

real*8 , allocatable :: dlong_deg(:), ddxyt(:,:,:), latmean(:), convRloc(:), &
                        dxrel(:,:), dyrel(:,:), dhrel(:,:), dtrel(:,:), &
			covar_scales(:,:,:), dd(:,:), exp_dd(:,:)


! constants used:
! -----------------
    pi19 = 3.141592653589793238 
    Rterre = 6371229.0
    radeg  = 180/pi19
    degrad = pi19/180
    convR = Rterre*degrad
    l_H = 5000
    n_gauss = 2
! -----------------


! ==========================================================
!          data-data correlation matrix  (F_dd)
! ==========================================================
 allocate (dlong_deg(nb_data), convRloc(nb_data), &
           latmean(nb_data), ddxyt(nb_data,nb_data,4))
 allocate (dxrel(nb_data,nb_data), dyrel(nb_data,nb_data), &
           dhrel(nb_data,nb_data), dtrel(nb_data,nb_data), &
           covar_scales(nb_data,nb_data,2))
 allocate (dd(nb_data,nb_data))

 ddxyt(1:nb_data,1:nb_data,1:4) = 1.0D12
 covar_scales(1:nb_data,1:nb_data,1:2) = 1.0

!  distances data-data :
 do j=1,nb_data  ! lower triangle
    latmean(1:j)  = (data_xyzt(1:j,2)+data_xyzt(j,2))*degrad*0.5    
    convRloc(1:j) = convR*cos(latmean(1:j))	!              15/10/2009
    dlong_deg(1:j) = data_xyzt(1:j,1) - data_xyzt(j,1)
    where (dlong_deg(1:j) .GE.  180) dlong_deg(1:j) = dlong_deg(1:j) - 360
    where (dlong_deg(1:j) .LE. -180) dlong_deg(1:j) = dlong_deg(1:j) + 360

    ddxyt(1:j,j,1) = dlong_deg(1:j)*convRloc(1:j)
    ddxyt(1:j,j,2) = (data_xyzt(1:j,2) - data_xyzt(j,2))*convR  
    ddxyt(1:j,j,3) = (data_xyzt(1:j,3) - data_xyzt(j,3))
    ddxyt(1:j,j,4) = (data_xyzt(1:j,4) - data_xyzt(j,4))
    do i = 1,j
       covar_scales(i,j,1) = min(abs(data_covms(i,1)),abs(data_covms(j,1)))
       covar_scales(i,j,2) = min(abs(data_covms(i,2)),abs(data_covms(j,2)))
    enddo
 enddo 

 
 do j=1,nb_data   ! upper triangle
    ddxyt(j+1:nb_data,j,1) = ddxyt(j,j+1:nb_data,1)
    ddxyt(j+1:nb_data,j,2) = ddxyt(j,j+1:nb_data,2)  
    ddxyt(j+1:nb_data,j,3) = ddxyt(j,j+1:nb_data,3) 
    ddxyt(j+1:nb_data,j,4) = ddxyt(j,j+1:nb_data,4)
    covar_scales(j+1:nb_data,j,1) = covar_scales(j,j+1:nb_data,1) 
    covar_scales(j+1:nb_data,j,2) = covar_scales(j,j+1:nb_data,2) 
 enddo 

 
! Scale 1 (uniform in space, x/y anisotropy)
   dxrel(1:nb_data,1:nb_data) = ddxyt(1:nb_data,1:nb_data,1)/covar_ls(1)
   dyrel(1:nb_data,1:nb_data) = ddxyt(1:nb_data,1:nb_data,2)/covar_ls(2)
   dtrel(1:nb_data,1:nb_data) = ddxyt(1:nb_data,1:nb_data,4)/covar_ls(3)
   dhrel(1:nb_data,1:nb_data) = ddxyt(1:nb_data,1:nb_data,3)/l_H
   dd(1:nb_data,1:nb_data)    = &
              config(1)*dxrel(1:nb_data,1:nb_data)*dxrel(1:nb_data,1:nb_data)&
            + config(2)*dyrel(1:nb_data,1:nb_data)*dyrel(1:nb_data,1:nb_data)&
            + config(3)*dhrel(1:nb_data,1:nb_data)*dhrel(1:nb_data,1:nb_data)&
            + config(4)*dtrel(1:nb_data,1:nb_data)*dtrel(1:nb_data,1:nb_data)
   F_dd(1:nb_data,1:nb_data,1) = exp(-dd(1:nb_data,1:nb_data)*0.5) 
   where(dd .GT. 100)  F_dd(:,:,1) = 0.0D0

! Scale 2 (location dependent, x/y anisotropy)
   dxrel(1:nb_data,1:nb_data) = & 
                ddxyt(1:nb_data,1:nb_data,1)/covar_scales(1:nb_data,1:nb_data,1)
   dyrel(1:nb_data,1:nb_data) = &
                ddxyt(1:nb_data,1:nb_data,2)/covar_scales(1:nb_data,1:nb_data,2)
   dtrel(1:nb_data,1:nb_data) = ddxyt(1:nb_data,1:nb_data,4)/covar_ms_t
   dhrel(1:nb_data,1:nb_data) = ddxyt(1:nb_data,1:nb_data,3)/l_H
   dd(1:nb_data,1:nb_data)    = &
               config(1)*dxrel(1:nb_data,1:nb_data)*dxrel(1:nb_data,1:nb_data)&
             + config(2)*dyrel(1:nb_data,1:nb_data)*dyrel(1:nb_data,1:nb_data)&
             + config(3)*dhrel(1:nb_data,1:nb_data)*dhrel(1:nb_data,1:nb_data)&
             + config(4)*dtrel(1:nb_data,1:nb_data)*dtrel(1:nb_data,1:nb_data)
   F_dd(1:nb_data,1:nb_data,2) = exp(-dd(1:nb_data,1:nb_data)*0.5)   
   where(dd .GT. 100)  F_dd(:,:,2) = 0.0D0

deallocate (dxrel, dyrel, dhrel, dtrel)
deallocate (ddxyt, latmean, convRloc, dlong_deg, dd, covar_scales)


! ==========================================================
!             covariance model-data (cmd)
! ==========================================================
allocate (ddxyt(nb_ana,nb_data,4),latmean(nb_ana), convRloc(nb_ana))
allocate (dxrel(nb_ana,nb_data), dyrel(nb_ana,nb_data), &
         dhrel(nb_ana,nb_data), dtrel(nb_ana,nb_data), &
	 covar_scales(nb_ana,nb_data,2))
allocate (dd(nb_ana,nb_data))


 ddxyt(1:nb_ana,1:nb_data,1:4) = 1.0D12
 covar_scales(1:nb_ana,1:nb_data,1:2) = 1.0


!  distances model-data:
!  ---------------------
do j=1,nb_data
      latmean(1:nb_ana) = (ana_xyz(1:nb_ana,2)+data_xyzt(j,2))*degrad*0.5
      convRloc(1:nb_ana) = convR*cos(latmean(1:nb_ana))
      ddxyt(1:nb_ana,j,1) = (ana_xyz(1:nb_ana,1) - &
                             data_xyzt(j,1))*convRloc(1:nb_ana)
      ddxyt(1:nb_ana,j,2) = (ana_xyz(1:nb_ana,2) - data_xyzt(j,2))*convR;
      ddxyt(1:nb_ana,j,3) =  ana_xyz(1:nb_ana,3) - data_xyzt(j,3);
      ddxyt(1:nb_ana,j,4) =  - data_xyzt(j,4);
 !     write(*,*) ana_covms(1,1:2), data_covms(j,1:2)
      do i = 1,nb_ana     
	 covar_scales(i,j,1) = min(abs(ana_covms(i,1)),abs(data_covms(j,1)))	 
         covar_scales(i,j,2) = min(abs(ana_covms(i,2)),abs(data_covms(j,2)))
     enddo
 !    write(*,*) covar_scales(1,j,1:2)
enddo


 ! Scale 1 (uniform in space, x/y anisotropy)
   dxrel(1:nb_ana,1:nb_data) = ddxyt(1:nb_ana,1:nb_data,1)/covar_ls(1)
   dyrel(1:nb_ana,1:nb_data) = ddxyt(1:nb_ana,1:nb_data,2)/covar_ls(2)
   dtrel(1:nb_ana,1:nb_data) = ddxyt(1:nb_ana,1:nb_data,4)/covar_ls(3)  
   dhrel(1:nb_ana,1:nb_data) = ddxyt(1:nb_ana,1:nb_data,3)/l_H
   dd(1:nb_ana,1:nb_data)    = config(1)*dxrel(1:nb_ana,1:nb_data)*dxrel(1:nb_ana,1:nb_data) &
                             + config(2)*dyrel(1:nb_ana,1:nb_data)*dyrel(1:nb_ana,1:nb_data) &
                             + config(3)*dhrel(1:nb_ana,1:nb_data)*dhrel(1:nb_ana,1:nb_data) &
                             + config(4)*dtrel(1:nb_ana,1:nb_data)*dtrel(1:nb_ana,1:nb_data)   
   F_md(1:nb_ana,1:nb_data,1) = exp(-dd(1:nb_ana,1:nb_data)*0.5)
   where(dd .GT. 100)  F_md(:,:,1) = 0.0D0
    
! Scale 2 (location dependent, isotropic)
   dxrel = ddxyt(1:nb_ana,1:nb_data,1)/covar_scales(1:nb_ana,1:nb_data,1)
   dyrel = ddxyt(1:nb_ana,1:nb_data,2)/covar_scales(1:nb_ana,1:nb_data,2)
   dtrel = ddxyt(1:nb_ana,1:nb_data,4)/covar_ms_t  
   dhrel = ddxyt(1:nb_ana,1:nb_data,3)/l_H
   dd    = config(1)*dxrel*dxrel + config(2)*dyrel*dyrel + &
           config(3)*dhrel*dhrel + config(4)*dtrel*dtrel
   F_md(1:nb_ana,1:nb_data,2) = exp(-dd(1:nb_ana,1:nb_data)*0.5)
   where(dd .GT. 100)  F_md(:,:,2) = 0.0D0

deallocate (dxrel, dyrel, dhrel, dtrel)
deallocate (ddxyt, latmean, convRloc, dd, covar_scales)

return
end


