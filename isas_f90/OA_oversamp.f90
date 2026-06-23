!==============================================================================
! 
 subroutine OA_oversamp (F_dd, nb_data, mx_gauss, var_weigh, w12, cov_max,  &
                         iok_data, nbok_data, iok_data2, nbok_data2)
!
!==============================================================================
!
!  Objet: Build list of data without oversampling (iok_data2)
!  -----
!  
!==============================================================================
!
! Inputs: 
integer*4    :: nb_data,mx_gauss, nbok_data, &
                iok_data(nbok_data)
real*8       :: var_weigh(3), & ! (relative weight of LS, MS, SG)
                cov_max, w12, F_dd(nb_data,nb_data,mx_gauss)  

! Outputs:
integer*4    :: nbok_data2, iok_data2(nbok_data)

! Locals:
integer*4   :: ii, nb_oversamp
real*8      :: xx

integer*4 , allocatable ::  listeij(:,:), mask(:) 
real*8 , allocatable   	::  xliste(:)

novmax = nbok_data*(nbok_data+1)/2
allocate (xliste(novmax), listeij(novmax,2))
 
 
!  locate oversampling data
!  ------------------------
listeij(1:nbok_data,1:2) = 0
nb_oversamp = 0
 
do j = 1, nbok_data 
   jj = iok_data(j)
   do i = 1, j-1 
       ii = iok_data(i)
       xx = (var_weigh(1)*F_dd(ii,jj,1) +  var_weigh(2)*F_dd(ii,jj,2))/w12  
       if (xx .GE. cov_max) then
          nb_oversamp = nb_oversamp + 1
          listeij(nb_oversamp,1) = i          
          listeij(nb_oversamp,2) = j
          xliste(nb_oversamp) = xx
       endif      
   enddo
enddo


! Create list of data to keep
!  ----------------------------
allocate (mask(nbok_data)) 

mask(1:nbok_data) = 1
do i = 1, nb_oversamp 
   mask(listeij(i,2)) = 0
enddo
nbok_data2 = sum(mask)

j = 0
do i = 1,nbok_data
   if (mask(i) .EQ. 1) then
      j = j + 1
      iok_data2(j) = iok_data(i)
   endif
enddo
nbok_data2 = j

!write(*,*) nbok_data2
!write(*,*) iok_data2(1:nbok_data2)


deallocate (xliste, listeij, mask)

return
end
