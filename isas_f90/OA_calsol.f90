
!==============================================================================
 
     subroutine  OA_calsol(C_md, im, C_dd, id, data_var, ana_var_pr, &
                              dino, nb_data, nb_ana, &
		              ana_fld, ana_var_ps, dres, cond_i, i_err)
 
!==============================================================================
!
!  Objet: Computes OA solution:
!  -----  
!   (X_est - X_0) = Koa (Y - Y_0)
!   Koa = C_md (C_dd + R)**-1
!   diag(C_est) =  ...
!
!   R = diag(data_var)
!   Y _ Y_O : dino
!   X_est - X_0 : ana_fld
!  diag(C_est)  : ana_var_ps
!
!   INPUT PARAMETERS: 
!   ----------------
!
!    C_md	: model-data cov matrix
!    im   	: leading dim of C_md
!    C_dd 	: data-data cov matrix
!    id 	: leading dim of C_dd
!    data_var	: data variance 
!    ana_var_pr : ap-priori variance on analyzed field
!    dino	: innovation value (Y - Y_0)
!    nb_data	: number of valid data
!    nb_ana	: number of analysis points
!       
!
!   OUTPUT PARAMETERS :
!   -----------------
!    ana_fld	: analyzed field 
!    ana_var_ps	: variance of analyzed field
!    dres	: residuals at data points
!    cond_i	: (C_dd + R) condition number
!    i_err 	: error indicator
!
!
!  Version: 4.0 f95 - mpi
!  -------
!  history
!  1.01  Creation d'apres programme fortran (Oct. 1990)	20/05/99  F.Gaillard
!  1.02  Revision filtrage moyenne echelle		22/12/99  F.Gaillard
!  1.03  Extraction calcul des distances		17/01/01  F.Gaillard
!  1.04  Ajout calcul des residus+homogeneisation 	23/11/01  F.Gaillard
!  V2.00 F. Gaillard, E. Autret - 			xx/06/2002 
!  V4.0	- F95-MPI 2					26/05/2006 F. Gaillard
! 
!==============================================================================

 implicit none


!  INPUTS:
!  ------
 integer*4  :: im, id, nb_data, nb_ana
 real*8     :: C_md(im,*), C_dd(id,*), data_var(*), ana_var_pr(*), dino(*)

!  OUTPUTS:
!  -------
 real*8     :: ana_fld(*), ana_var_ps(*), &
               dres(*),    &
	       cond_i	
 integer*4  :: i_err

!  LOCALS:
! -------
 integer*4  ::i, j
 
 character*1  ::  tbla, tblb, side, uplo	!  LINPACK/BLAS
 real*8       ::  abla, bbla, Cdd_norm, cov_min
 external     ::  dpotrs, dpocon, dgemm
 integer*4    ::  incx, incy
 
real*4 tps_cal, tps_beg, tps_end 

		
!  DYNAMICALLY ALLOCATED ARRAYS

 real*8 ,   allocatable :: WKSP0(:,:), WKSP1(:), WKSP2(:,:), C_dd_tot(:,:)
 integer*4, allocatable :: IWORK(:)
 
 !  Program starts

! test !   call cpu_time(tps_beg)
 
 cond_i = 0
 cov_min = 0.000001
  
!  Builds matrix to invert:
!  ----------------------
 allocate (C_dd_tot(nb_data, nb_data))
 C_dd_tot(1:nb_data,1:nb_data) = C_dd(1:nb_data,1:nb_data)
 where (C_dd_tot(1:nb_data,1:nb_data) .LT. cov_min) &
            C_dd_tot =  0.0000000000000D0
   
 do i = 1,  nb_data
    C_dd_tot(i,i) = C_dd_tot(i,i) + data_var(i)
 enddo

 
 
! Decomposition de DD en triangulaire + diagonale (Cholesky)
! ----------------------------------------------------------
 i_err  = 0
 call dpotrf('U', nb_data, C_dd_tot, nb_data, i_err)
 if (i_err .ne. 0) then
    return
 endif
 
 allocate (WKSP1(3*nb_data), IWORK(nb_data))
 Cdd_norm = 1
 
 call dpocon('U', nb_data, C_dd_tot, nb_data, Cdd_norm, cond_i, WKSP1, IWORK, i_err)
 deallocate (WKSP1, IWORK)
 
 
! Resolution
! ----------
 allocate (WKSP1(nb_data))
 WKSP1(1:nb_data) = dino(1:nb_data)
 
 call dpotrs('U', nb_data, 1, C_dd_tot, nb_data, WKSP1, nb_data, i_err)
 if (i_err .ne. 0) return


! =================================
!     Calcul des residus
! =================================
  do i = 1, nb_data
     dres(i) = data_var(i)*WKSP1(i)
  end do


! =======================================================
!     Calcul de l'estimateur: (X_est - X_0) = Koa (Y - Y_0)
!                              Koa = C_md (C_dd + R)**-1
! =======================================================
   abla = 1.d0			! parametres pour la blas
   bbla = 1.d0
   tbla = 'n'
   incx = 1
   incy = 1

   ana_fld(1: nb_ana) = 0

   call dgemv (tbla, nb_ana, nb_data, abla, C_md, im, WKSP1, incx, &
                     bbla, ana_fld, incy)

 ! test !   call cpu_time(tps_end)
 ! test !  tps_cal = tps_end - tps_beg
 ! test !  write(*,*) 'estim' , tps_cal
 ! test !  call cpu_time(tps_beg)
  
!  ========================================	
!            Computes error:
!  ========================================
!		     
!  calcule DD-1.*CmdT:
!  -----------------
   allocate (WKSP0(nb_data, nb_ana), WKSP2(nb_data, nb_ana))

   WKSP0 = transpose(C_md(1:nb_ana,1:nb_data))
   WKSP2 = WKSP0
   uplo = 'U'
   call dpotrs(uplo, nb_data, nb_ana ,C_dd_tot, nb_data, WKSP0, nb_data, i_err)
   
 
  
!  Computes diagonal of error covariance matrix
   do i = 1, nb_ana
      ana_var_ps(i) = 0
      do j = 1, nb_data
         ana_var_ps(i) = ana_var_ps(i) + WKSP2(j,i)*WKSP0(j,i)
      enddo
      ana_var_ps(i) = ana_var_pr(i) - ana_var_ps(i)
   enddo

 deallocate (WKSP0, WKSP1, WKSP2, C_dd_tot)
 
 ! test !   call cpu_time(tps_end)
 ! test !   tps_cal = tps_end - tps_beg
 ! test !   write(*,*) 'error' , tps_cal
 
return
end
