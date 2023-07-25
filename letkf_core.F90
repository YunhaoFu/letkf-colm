!--F90

!#########################################################################
!# File Name: letkf_core.F90
!# Created on : 2022-12-03 01:12:05
!# Author: YunhaoFu
!# Last Modified: 2023-04-13 09:06:36
!# Description: works for fortran90 with pre-compile
!#########################################################################

MODULE letkf_core
    USE precision
    IMPLICIT NONE

    CONTAINS

          !================================================================================
          !> The core LETKF equation
          !--------------------------------------------------------------------------------
          SUBROUTINE letkf_core_solve(nbv, nobs, hdxb, rdiag, rloc, dep, infl, trans)
            INTEGER, INTENT(in)   :: nbv             !< ens_size
            INTEGER, INTENT(in)   :: nobs            !< number of observations
            REAL(r8), INTENT(in)  :: hdxb(nbv,nobs)  !< ensemble perturbations in obs space
            REAL(r8), INTENT(in)  :: rdiag(nobs)     !< observation error variance
            REAL(r8), INTENT(in)  :: rloc(nobs)      !< observation localization weights
            REAL(r8), INTENT(in)  :: dep(nobs)       !< observation departures
            REAL(r8), INTENT(in)  :: infl            !< covariance inflation
            REAL(r8), INTENT(out) :: trans(nbv, nbv) !< ensemble transformation matrix

            ! temporary intermediate values
            REAL(r8)              :: hdxb_rinv(nbv,nobs)
            REAL(r8)              :: work1(nbv,nbv)
            REAL(r8)              :: work2(nbv,nobs)
            REAL(r8)              :: work3(nbv)
            REAL(r8)              :: pa(nbv,nbv)
            REAL(r8)              :: eival(nbv)
            REAL(r8)              :: eivec(nbv,nbv)
            REAL(r8), ALLOCATABLE :: evwork(:)
            INTEGER               :: i, j, err
            REAL(r8)              :: r

            ALLOCATE(evwork((64+2)*nbv))

            ! hdxb Rinv
            DO j = 1, nobs
               hdxb_rinv(:,j) = hdxb(:,j) / rdiag(j) * rloc(j)
            END DO

            ! hdxb^T Rinv hdxb
            CALL dgemm('n','t', nbv, nbv, nobs, &
                 1.0_r8, hdxb_rinv, nbv, hdxb, nbv, 0.0_r8, work1, nbv)


            ! hdxb^T Rinv hdxb + (k-1) I / rho (covariance inflation)
            r = REAL(nbv-1,r8)*1.0_r8/infl
            DO i=1,nbv
               work1(i,i) = work1(i,i) + r
            END DO

            ! eigenvalues and eigenvectors of above
            CALL dsyev('V','U', nbv, work1, nbv, eival, evwork, SIZE(evwork), err)
            eivec = work1

            ! Pa = [hdxb^T Rinv hdxb + (m-1) I] inv
            DO i=1,nbv
               work1(:,i) = eivec(:,i) / eival(i)
            END DO
            CALL dgemm('n','t',nbv,nbv,nbv,1.0_r8, work1, nbv, eivec,&
                 nbv, 0.0_r8, pa, nbv)

            ! Pa hdxb_rinv^T
            CALL dgemm('n', 'n', nbv, nobs, nbv, 1.0_r8, pa, nbv, hdxb_rinv,&
                 nbv, 0.0_r8, work2, nbv)

            ! Pa hdxb_rinv^T dep
            work3 = 0
            DO j=1,nobs
               work3 = work3 + work2(:,j) * dep(j)
            END DO

            ! T = sqrt[(m-1)Pa]
            DO j = 1, nbv
               r = SQRT(REAL(nbv-1,r8) / eival(j))
               work1(:,j) = eivec(:,j) * r
            END DO
            CALL dgemm('n', 't', nbv, nbv, nbv, 1.0_r8, work1, nbv, eivec, &
                 nbv, 0.0_r8, trans, nbv)

            ! Relaxation (RTPP or RTPS?)
            ! TODO

            ! T + Pa hdxb_rinv^T dep
            DO j=1,nbv
               trans(:,j) = trans(:,j) + work3
            END DO

            ! adaptive inflation
            ! TODO
            DEALLOCATE(evwork)

          END SUBROUTINE letkf_core_solve
          !================================================================================

END MODULE letkf_core
