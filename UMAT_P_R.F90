!                                                                      C
!     UMAT_PCLI_R.for     Plasticity Classical Theory                  C
!                                                                      C
!     LOCAL ARRAYS                                                     C
!                                                                      C
!     EELAS - ELASTIC STRAINS STATEV(1..NTENS)                         C
!     EPLAS - PLASTIC STRAINS STATEV(NTENS+1...2*NTENS)                C
!     FLOW  - PLASTIC FLOW DIRECTION                                   C
!     HARD  - HARDENING MODULUS                                        C
!                                                                      C
!     PROPS(1) - E                                                     C
!     PROPS(2) - NU                                                    C
!     PROPS(5..) - SYIELD AN HARDENING DATA                            C
!     CALLS KUHARD FOR CURVE OF YIELD STRESS VS. PLASTIC STRAINS       C
!     EQPLAS - EQUIVALENT PLASTIC STRAIN STATEV(2*NTENS+1)             C
!                                                                      C

      SUBROUTINE UMAT(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL,&
      DDSDDT, DRPLDE, DRPLDT, STRAN, DSTRAN, TIME, DTIME, TEMP,&
      DTEMP, PREDEF, DPRED, CMNAME, NDI, NSHR, NTENS, NSTATV,&
      PROPS, NPROPS, COORDS, DROT, PNEWDT, CELENT, DFGRD0,&
      DFGRD1, NOEL, NPT, LAYER, KSPT, KSTEP, KINC)

      IMPLICIT REAL*8(A-H,O-Z)

      CHARACTER*80 CMNAME
 
      DIMENSION STRESS(NTENS), STATEV(NSTATV),DDSDDE(NTENS, NTENS),&
      DDSDDT(NTENS),STRAN(NTENS),DSTRAN(NTENS),TIME(2),PROPS(NPROPS) 
     

      DIMENSION EELAS(NTENS), EPLAS(NTENS),DS(NTENS),DSTRESS(NTENS)


      DIMENSION AUX1(1,NTENS),AUX2(NTENS,NTENS),AUX3(NTENS,NTENS),&
      AUX4(NTENS,NTENS),AUX5(NTENS,NTENS),AUX6(NTENS,1),AUX7(1,NTENS),&
      AUX8(1,1),AUX9(1,NTENS),AUX10(1,1),AUX11(NTENS,NTENS),&
      AUX12(NTENS,1),AUX13(1,NTENS),AUX14(1,1),AUX15(1,NTENS),&
      AUX16(NTENS,NTENS),AUX17(NTENS,NTENS),STRESST(1,NTENS),&
      P(NTENS,NTENS),SINVAR(1,1),BI(NTENS,NTENS),STRESSUPD(NTENS,1),&
      SDEV(NTENS,1),EM(NTENS,NTENS),B(NTENS,NTENS),Q(NTENS,NTENS),&
      QT(NTENS,NTENS),DP(NTENS,NTENS),DC(NTENS,NTENS),DEL(NTENS,NTENS),&
      BT(NTENS,NTENS),DIAG(NTENS,NTENS),GDIA(NTENS,NTENS)

      PARAMETER (ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0, SIX=6.0D0,&
                 ENUMAX=0.4999D0, NEWTON=10, TOLER=1.0D-7,MAXITER=30,&
                 FOUR=4.D0)

!     Inititalize arrays

      CALL KCLEAR(AUX1,1,NTENS)
      CALL KCLEAR(AUX2,NTENS,NTENS)
      CALL KCLEAR(AUX3,NTENS,NTENS)
      CALL KCLEAR(AUX4,NTENS,NTENS)
      CALL KCLEAR(AUX5,NTENS,NTENS)
      CALL KCLEAR(AUX6,NTENS,1)
      CALL KCLEAR(AUX7,1,NTENS)
      CALL KCLEAR(AUX8,1,1)
      CALL KCLEAR(AUX9,1,NTENS)
      CALL KCLEAR(AUX10,1,1)
      CALL KCLEAR(AUX11,NTENS,NTENS)
      CALL KCLEAR(AUX12,NTENS,1)
      CALL KCLEAR(AUX13,1,NTENS)
      CALL KCLEAR(AUX14,1,1)
      CALL KCLEAR(AUX15,1,NTENS)
      CALL KCLEAR(AUX16,NTENS,NTENS)
      CALL KCLEAR(AUX17,NTENS,NTENS)

      CALL KCLEAR(P,NTENS,NTENS)
      CALL KCLEAR(Q,NTENS,NTENS)
      CALL KCLEAR(DP,NTENS,NTENS)
      CALL KCLEAR(DC,NTENS,NTENS)
      CALL KCLEAR(QT,NTENS,NTENS)
      CALL KCLEAR(DEL,NTENS,NTENS)

      CALL KCLEAR(STRESST,1,NTENS)
      CALL KCLEAR(SINVAR,1,1)

      CALL KCLEAR(STRESSUPD,NTENS,1)
      CALL KCLEAR(AUX17,NTENS,NTENS)
      CALL KCLEAR(B,NTENS,NTENS)
      CALL KCLEAR(BT,NTENS,NTENS)
 
!     Recover equivalent plastic strain, elastic strains, and plastic
!     strains. Also initialize user definde data sets.

      DO K1=1, NTENS
        EELAS(K1)=STATEV(K1)
        EPLAS(K1)=STATEV(K1+NTENS)
      END DO
      EQPLAS=STATEV(1+2*NTENS)

!     Elastic properties

      EMOD=PROPS(1)
      ENU=MIN(PROPS(2),ENUMAX)
      EBULK3=EMOD/(ONE-TWO*ENU)
      EG2=EMOD/(ONE+ENU)
      EG=EG2/TWO
      EG3=THREE*EG
      ELAM=(EBULK3-EG2)/THREE
      CBETA2=EG2

!     Elastic stiffness

      DO K1=1, 3
        DO K2=1, 3
          DDSDDE(K2, K1)=ELAM
        END DO
        DDSDDE(K1, K1)=EG2+ELAM
      END DO
      DO K1=4, 4
        DDSDDE(K1, K1)=EG
      END DO

!     Form proyector and diagonal decomposition matrices

      CALL KPROYECTOR(P)

!     Calculate predictor stress and elastic strains

      CALL KMAVEC(DDSDDE,NTENS,NTENS,DSTRAN,DS)
      CALL KUPDVEC(STRESS,NTENS,DS)
      CALL KUPDVEC(EELAS,NTENS,DSTRAN)

!     Calculate equivalent mises stress

      CALL KMTRAN(STRESS,NTENS,1,STRESST)
      CALL KMMULT(P,NTENS,NTENS,STRESS,NTENS,1,SDEV)
      CALL KMMULT(STRESST,1,NTENS,SDEV,NTENS,1,SINVAR)
      FBAR=DSQRT(SINVAR(1,1))

!     Get the yield stress from the specifid hardening function.

      CALL KUHARD(SYIEL0,EHARD,EQPLAS,2,PROPS(3))

!     Determine if actively yielding

      SYIELD=SYIEL0
      IF(FBAR.GT.(ONE+TOLER)*SYIEL0) THEN

!       Actively yielding-Perform local Newton iterations
!       to find consistncy parameter and equivalent plastic
!       strain

!       Starts iterations

        ITER=1
        GAM_PAR=ZERO
        IFLAG=0
        DO
          CALL KSPECTRAL(Q,DP,DC,DIAG,GDIA,GAM_PAR,EMOD,ENU)
          CALL KMTRAN(Q,NTENS,NTENS,QT)
          CALL KMMULT(STRESST,1,NTENS,Q,NTENS,NTENS,AUX1)
          CALL KMMULT(QT,NTENS,NTENS,P,NTENS,NTENS,AUX2)
          CALL KMMULT(AUX2,NTENS,NTENS,Q,NTENS,NTENS,AUX3)
          CALL KMMULT(AUX3,NTENS,NTENS,DIAG,NTENS,NTENS,AUX4)
          CALL KMMULT(AUX4,NTENS,NTENS,QT,NTENS,NTENS,AUX5)
          CALL KMMULT(AUX5,NTENS,NTENS,STRESS,NTENS,1,AUX6)
          CALL KMMULT(AUX1,1,NTENS,DIAG,NTENS,NTENS,AUX7)
          CALL KMMULT(AUX7,1,NTENS,AUX6,NTENS,1,AUX8)
          CALL KMMULT(AUX1,1,NTENS,GDIA,NTENS,NTENS,AUX9)
          CALL KMMULT(AUX9,1,NTENS,AUX6,NTENS,1,AUX10)
          FBAR=DSQRT(AUX8(1,1))
          TETA2=ONE-(TWO/THREE)*EHARD*GAM_PAR
          FJAC=TETA2*AUX10(1,1)/FBAR-(TWO/THREE)*EHARD*FBAR
          FGAM=FBAR-SYIELD

!         Updates

          GAM_PAR=GAM_PAR-FGAM/FJAC
          EQPLAS1=EQPLAS+DSQRT(TWO/THREE)*GAM_PAR*FBAR
          CALL KUHARD(SYIELD,EHARD,EQPLAS1,2,PROPS(3))

          IF(ABS(FGAM/FJAC).LT.TOLER) THEN
            IFLAG=0
            GOTO 801
          ELSE
            IF(ITER.GT.MAXITER) THEN
              IFLAG=1
              GOTO 802
            END IF
          END IF

          ITER=ITER+1
        END DO

  801   CONTINUE

!       Local Newton algorithm converged
!       Update stresses, elastic and plastic strains, equivalent plastic
!       strains

        CALL KSPECTRAL(Q,DP,DC,DIAG,GDIA,GAM_PAR,EMOD,ENU)
        CALL KMMULT(Q,NTENS,NTENS,DIAG,NTENS,NTENS,AUX17)
        CALL KMMULT(AUX17,NTENS,NTENS,QT,NTENS,NTENS,B)
        CALL KMMULT(B,NTENS,NTENS,STRESS,NTENS,1,STRESSUPD)

        CALL KCLEAR(STRESS,NTENS,1)
        DO K1=1,NTENS
          STRESS(K1)=STRESSUPD(K1,1)
        END DO

        CALL KCLEAR(STRESST,1,NTENS)
        CALL KCLEAR(SDEV,NTENS,1)
        CALL KMTRAN(STRESS,NTENS,1,STRESST)
        CALL KMMULT(P,NTENS,NTENS,STRESS,NTENS,1,SDEV)
        CALL KMMULT(STRESST,1,NTENS,SDEV,NTENS,1,SINVAR)
        FBAR=DSQRT(SINVAR(1,1))

        DO K1=1,NTENS
          EPLAS(K1)=EPLAS(K1)+GAM_PAR*SDEV(K1,1)
          EELAS(K1)=EELAS(K1)-EPLAS(K1)
        END DO

        EQPLAS=EQPLAS1

!       Formulate the consistent material Jacobian (tangent)

        CALL KCLEAR(EM,NTENS,NTENS)
        CALL KMMULT(B,NTENS,NTENS,DDSDDE,NTENS,NTENS,EM)
        CALL KMMULT(EM,NTENS,NTENS,P,NTENS,NTENS,AUX11)
        CALL KMMULT(AUX11,NTENS,NTENS,STRESS,NTENS,1,AUX12)
        CALL KMMULT(STRESST,1,NTENS,P,NTENS,NTENS,AUX13)
        CALL KMMULT(AUX13,1,NTENS,AUX12,NTENS,1,AUX14)
        CALL KMTRAN(AUX12,NTENS,1,AUX15)
        CALL KMMULT(AUX12,NTENS,1,AUX15,1,NTENS,AUX16)
        SCALAR1=ONE/AUX14(1,1)
        CALL KSMULT(AUX16,NTENS,NTENS,SCALAR1)
        TETA2=ONE-(TWO/THREE)*EHARD*GAM_PAR
        CBETA=(TWO/THREE/TETA2/AUX14(1,1))*FBAR*FBAR*EHARD
        SCALAR2=ONE/(ONE+CBETA)
        CALL KSMULT(AUX16,NTENS,NTENS,SCALAR2)
        CALL KCLEAR(DDSDDE,NTENS,NTENS)
        CALL KMATSUB(EM,NTENS,NTENS,AUX16,DDSDDE,0)

      END IF

!     Store elastic strains, (equivalent) plastic strains
!     in state variable array

      DO K1=1,NTENS
        STATEV(      K1)=EELAS(K1)
        STATEV(NTENS+K1)=EPLAS(K1)
      END DO
      STATEV(2*NTENS+1)=EQPLAS
      STATEV(2*NTENS+2)=DSQRT(THREE/TWO)*FBAR

  802 IF (IFLAG.EQ.1) THEN
         WRITE(*,*)
         WRITE(*,*) 'LOCAL PLASTICITY ALGORITHM DID NOT CONVREGED'
         WRITE(*,*) 'AT GAUSS POINT=',NPT, 'ELEMENT=',NOEL
         WRITE(*,*) 'AFTER=',ITER,' ITERATIONS'
         WRITE(*,*) 'LAST CORRECTION=',FGAM/FJAC
         CALL XIT
      END IF

      END SUBROUTINE UMAT

!                                                                      C
!     SUBROUTINE UHARD                                                 C
!                                                                      C

      SUBROUTINE KUHARD(SYIELD,EHARD,EQPLAS,NVALUE,TABLE)

      IMPLICIT REAL*8(A-H,O-Z)

      DIMENSION TABLE(2,NVALUE)

      PARAMETER(ZERO=0.D0,TWO=2.D0,THREE=3.D0)

      SYIEL0=TABLE(1,1)

!     Compute hardening modulus

      EHARD=(TABLE(1,2)-TABLE(1,1))/TABLE(2,2)

!     Compute yield stress corresponding to EQPLAS

      SYIELD=DSQRT(TWO/THREE)*(SYIEL0+EHARD*EQPLAS)

      RETURN

      END

!                                                                      C
!     SUBROUTINE SPECTRAL                                              C
!                                                                      C

      SUBROUTINE KSPECTRAL(Q,DP,DC,DIAG,GDIA,GAM_PAR,E,ENU)

      IMPLICIT REAL*8(A-H,O-Z)

      PARAMETER(ZERO=0.D0,ONE=1.D0,TWO=2.D0,THREE=3.D0,FOUR=4.D0,&
                SIX=6.D0)

      DIMENSION Q(4,4),DP(4,4),DC(4,4),DIAG(4,4),GDIA(4,4)

      CALL KCLEAR(Q,4,4)
      CALL KCLEAR(DP,4,4)
      CALL KCLEAR(DC,4,4)
      CALL KCLEAR(DIAG,4,4)
      CALL KCLEAR(GDIA,4,4)

      EG2=E/(ONE+ENU)
      EG=EG2/TWO
      ELAM=EG2*ENU/(ONE-TWO*ENU)
      CBETA2=EG2
      CALFA2=EG2

      Q(1,1)=ZERO
      Q(1,2)=TWO/DSQRT(SIX)
      Q(1,3)=ONE/DSQRT(THREE)
      Q(2,1)=-DSQRT(TWO)/TWO
      Q(2,2)=-ONE/DSQRT(SIX)
      Q(2,3)=ONE/DSQRT(THREE)
      Q(3,1)=DSQRT(TWO)/TWO
      Q(3,2)=-ONE/DSQRT(SIX)
      Q(3,3)=ONE/DSQRT(THREE)
      Q(4,4)=ONE

      DP(1,1)=ONE
      DP(2,2)=ONE
      DP(3,3)=ZERO
      DP(4,4)=TWO

      DC(1,1)=EG2
      DC(2,2)=EG2
      DC(3,3)=THREE*ELAM+EG2
      DC(4,4)=EG

      DIAG(1,1)=ONE/(ONE+CBETA2*GAM_PAR)
      DIAG(2,2)=ONE/(ONE+CBETA2*GAM_PAR)
      DIAG(3,3)=ONE
      DIAG(4,4)=ONE/(ONE+CBETA2*GAM_PAR)
      GDIA(1,1)=-(CBETA2/((ONE+CBETA2*GAM_PAR)**2))
      GDIA(2,2)=-(CBETA2/((ONE+CBETA2*GAM_PAR)**2))
      GDIA(3,3)=ZERO
      GDIA(4,4)=-(CBETA2/((ONE+CBETA2*GAM_PAR)**2))

      RETURN

      END

!                                                                      C
!     SUBROUTINE PROYECTOR                                             C
!                                                                      C
      SUBROUTINE KPROYECTOR(P)

      IMPLICIT REAL*8(A-H,O-Z)

      PARAMETER(ZERO=0.D0,ONE=1.D0,TWO=2.D0,THREE=3.D0)

      DIMENSION P(4,4)

      CALL KCLEAR(P,4,4)

      P(1,1)=TWO/THREE
      P(1,2)=-ONE/THREE
      P(1,3)=-ONE/THREE
      P(2,1)=-ONE/THREE
      P(2,2)=TWO/THREE
      P(2,3)=-ONE/THREE
      P(3,1)=-ONE/THREE
      P(3,2)=-ONE/THREE
      P(3,3)=TWO/THREE
      P(4,4)=TWO

      RETURN

      END


!                                                                      C
!             M A T R I X   H A N D L I N G                            C
!-------------U T I L I T I E S   B L O C K--------------              C
!                                                                      C
!                                                                      C
!      SUBROUTINE KCLEAR(A,N,M)                                        C
!      Clear a real matrix                                             C
!                                                                      C
      SUBROUTINE KCLEAR(A,N,M)

      IMPLICIT REAL*8(A-H,O-Z)

      PARAMETER(ZERO=0.0D0)
      DIMENSION A(N,M)

      DO I=1,N
        DO J=1,M
          A(I,J)=ZERO
        END DO
      END DO

      RETURN

      END
!                                                                      C
!      SUBROUTINE KMMULT(A,NRA,NCA,B,NRB,NCB,C)                        C
!      Real matrix product                                             C
!                                                                      C
      SUBROUTINE KMMULT(A,NRA,NCA,B,NRB,NCB,C)

      IMPLICIT REAL*8 (A-H,O-Z)

      PARAMETER(ZERO=0.D0)
      DIMENSION A(NRA,NCA),B(NRB,NCB),C(NRA,NCB)

      CALL KCLEAR(C,NRA,NCB)
      DUM=ZERO
      DO I=1,NRA
        DO J=1,NCB
         DO K=1,NCA
           DUM=DUM+A(I,K)*B(K,J)
          END DO
          C(I,J)=DUM
          DUM=ZERO
        END DO
      END DO

      RETURN

      END
!                                                                      C
!      SUBROUTINE KSMULT(A,NR,NC,S)                                    C
!      Matrix times a scalar.                                          C
!                                                                      C
      SUBROUTINE KSMULT(A,NR,NC,S)

      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION A(NR,NC)

      DO I=1,NR
        DO J=1,NC
          DUM=A(I,J)
          A(I,J)=S*DUM
          DUM=0.D0
        END DO  
      END DO

      RETURN

      END
!                                                                      C
!      SUBROUTINE KUPDMAT(A,NR,NC,B)                                   C
!      Updates an existing matrix with an incremental matrix.          C
!                                                                      C
      SUBROUTINE KUPDMAT(A,NR,NC,B)

      IMPLICIT REAL*8 (A-H,O-Z)

      PARAMETER(ZERO=0.D0)

      DIMENSION A(NR,NC),B(NR,NC)

      DO I=1,NR
        DO J=1,NC
          DUM=A(I,J)
          A(I,J)=ZERO
          A(I,J)=DUM+B(I,J)
          DUM=ZERO
        END DO
      END DO

      RETURN

      END
!                                                                      C
!      SUBROUTINE KMTRAN(A,NRA,NCA,B)                                  !      
!      Matrix transpose                                                C
!                                                                      C
      SUBROUTINE KMTRAN(A,NRA,NCA,B)

      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION A(NRA,NCA),B(NCA,NRA)

      CALL KCLEAR(B,NCA,NRA)
      DO I=1,NRA
       DO J=1,NCA
         B(J,I)=A(I,J)
        END DO
      END DO

      RETURN

      END
!                                                                      C
!      SUBROUTINE KMAVEC(A,NRA,NCA,B,C)                                C
!      Real matrix times vector                                        C
!                                                                      C
      SUBROUTINE KMAVEC(A,NRA,NCA,B,C)

      IMPLICIT REAL*8 (A-H,O-Z)

      PARAMETER(ZERO=0.D0)
      DIMENSION A(NRA,NCA),B(NCA),C(NRA)

      CALL KCLEARV(C,NRA)

      DO K1=1,NRA
        DO K2=1,NCA
          C(K1)=C(K1)+A(K1,K2)*B(K2)	    
        END DO
      END DO     

      RETURN

      END
!                                                                      C
!      SUBROUTINE KCLEARV(A,N)                                         C
!      Clear a real vector                                             C
!                                                                      C
      SUBROUTINE KCLEARV(A,N)

      IMPLICIT REAL*8(A-H,O-Z)

      PARAMETER(ZERO=0.0D0)

      DIMENSION A(N)

      DO I=1,N
        A(I)=ZERO
      END DO

      RETURN

      END
!                                                                      C
!      SUBROUTINE KUPDVEC(A,NR,B)                                      C
!      Updates an existing vector with an incremental vector.          C
!                                                                      C
      SUBROUTINE KUPDVEC(A,NR,B)

      IMPLICIT REAL*8 (A-H,O-Z)

      PARAMETER(ZERO=0.D0)

      DIMENSION A(NR),B(NR)

      DO I=1,NR
        DUM=A(I)
        A(I)=ZERO
        A(I)=DUM+B(I)
        DUM=ZERO
      END DO

      RETURN

      END

!                                                                      C
!      SUBROUTINE KVECSUB(A,NRA,B,NRB,C)                               C
!      Substracts one column vector from another column vector         C
!      IFLAG=0 for substraction                                        C
!      IFLAG=1 for addition                                            C
!                                                                      C
      SUBROUTINE KVECSUB(A,NRA,B,NRB,C,IFLAG)

      IMPLICIT REAL*8(A-H,O-Z)

      PARAMETER (ONE=1.0D0, ONENEG=-1.0D0)

      DIMENSION A(NRA,1),B(NRB,1),C(NRB,1)

      SCALAR=ONENEG

      IF (IFLAG.EQ.1) SCALAR=ONE

      DO I=1,NRA
        C(I,1)=A(I,1)+B(I,1)*SCALAR
      END DO

      RETURN

      END

!                                                                      C
!      SUBROUTINE KMATSUB(A,NRA,NCA,B,C,IFLAG)                         C
!      Substracts one rectangular matrix from another rectangular      C
!      matrix                                                          C
!      IFLAG=0 for substraction                                        C
!      IFLAG=1 for addition                                            C
!                                                                      C
      SUBROUTINE KMATSUB(A,NRA,NCA,B,C,IFLAG)

      IMPLICIT REAL*8(A-H,O-Z)

      PARAMETER (ONE=1.0D0, ONENEG=-1.0D0)

      DIMENSION A(NRA,NCA),B(NRA,NCA),C(NRA,NCA)

      CALL KCLEAR(C,NRA,NCA)

      SCALAR=ONENEG

      IF (IFLAG.EQ.1) SCALAR=ONE

      DO I=1,NRA
        DO J=1,NCA
          C(I,J)=A(I,J)+B(I,J)*SCALAR
        END DO
      END DO

      RETURN

      END


!                                                                      C
!     SUBROUTINE IDENTITY                                              C
!     CREATES AN IDENTITY MATRIX OF DIMENSIONS NDIM,NDIM               C
!                                                                      C

      SUBROUTINE KIDENTITY(DEL,NDIM)

      IMPLICIT REAL*8(A-H,O-Z)

      PARAMETER(ONE=1.D0)

      DIMENSION DEL(NDIM,NDIM)

      CALL KCLEAR(DEL,NDIM,NDIM)

      DO K1=1,NDIM
        DEL(K1,K1)=ONE
      END DO

      RETURN

      END

!                                                                      C
!   SUBROUTINE KINVERSE                                                 C
!                                                                      C
!   IVEERSE OF A MATRIX USING LU DECOMPOSITION                         C
!   TAKEN FROM NUMERICAL RECIPES By Press et al                        C
!                                                                      C
!   A   Matrix to be inverted.                                         C
!   Y   Inverse of A                                                   C
!   N   Dimension                                                      C
!                                                                      C
!                                                                      C
!                                                                      C

      SUBROUTINE KINVERSE(A,Y,NP,N)

      IMPLICIT REAL*8(A-H,O-Z)

      PARAMETER (ZERO=0.D0,ONE=1.D0)

      DIMENSION A(NP,NP),Y(NP,NP),INDX(NP),AUX(NP,NP)

      CALL KCLEAR(AUX,NP,NP)
      CALL KCOPYMAT(A,AUX,N)

      DO I=1,N
        DO J=1,N
          Y(I,J)=ZERO
        END DO
        Y(I,I)=ONE
      END DO
      CALL KLUDCMP(AUX,N,NP,INDX,D)
      DO J=1,N
        CALL KLUBKSB(AUX,N,NP,INDX,Y(1,J))
      END DO

      RETURN

      END

!                                                                      C
!   SUBROUTINE KLUDCMP                                                 C
!                                                                      C
!   LU MATRIX DECOMPOSITION                                            C
!   TAKEN FROM NUMERICAL RECIPES By Press et al                        C
!                                                                      C
!                                                                      C
!                                                                      C
      SUBROUTINE KLUDCMP(A,N,NP,INDX,D)

      IMPLICIT REAL*8(A-H,O-Z)

      PARAMETER(NMAX=500,TINY=1.0E-20,ZERO=0.D0,ONE=1.D0)

      DIMENSION INDX(N),A(NP,NP),VV(NMAX)

      D=ONE
      DO I=1,N
        AAMAX=ZERO
        DO J=1,N
          IF(ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
        END DO
        IF(AAMAX.EQ.0.) write(*,*) 'SINGULAR MATRIX IN LUDCMP'
        VV(I)=ONE/AAMAX
      END DO

      DO J=1,N
        DO I=1,J-1
          SUM=A(I,J)
          DO K=1,I-1
            SUM=SUM-A(I,K)*A(K,J)
          END DO
          A(I,J)=SUM
        END DO
        AAMAX=ZERO
        DO I=J,N
          SUM=A(I,J)
          DO K=1,J-1
            SUM=SUM-A(I,K)*A(K,J)
          END DO
          A(I,J)=SUM
          DUM=VV(I)*ABS(SUM)
          IF(DUM.GE.AAMAX) THEN
             IMAX=I
             AAMAX=DUM
          END IF
        END DO
        IF(J.NE.IMAX) THEN
           DO K=1,N
             DUM=A(IMAX,K)
             A(IMAX,K)=A(J,K)
             A(J,K)=DUM
           END DO
           D=-D
           VV(IMAX)=-VV(J)
        END IF
        INDX(J)=IMAX
        IF(A(J,J).EQ.0.) A(J,J)=TINY
        IF(J.NE.N) THEN
           DUM=ONE/A(J,J)
           DO I=J+1,N
             A(I,J)=A(I,J)*DUM
           END DO
        END IF
      END DO

      RETURN

      END

!                                                                      C
!   SUBROUTINE KLUBKSB                                                 C
!                                                                      C
!   FORWARD SUBSTITUTION                                               C
!   TAKEN FROM NUMERICAL RECIPES By Press et al                        C
!                                                                      C
!                                                                      C
!                                                                      C
      SUBROUTINE KLUBKSB(A,N,NP,INDX,B)

      IMPLICIT REAL*8(A-H,O-Z)

      PARAMETER (ZERO=0.D0)

      DIMENSION INDX(N),A(NP,NP),B(NP)

      II=0
      DO I=1,N
        LL=INDX(I)
        SUM=B(LL)
        B(LL)=B(I)
        IF(II.NE.0) THEN
           DO J=II,I-1
             SUM=SUM-A(I,J)*B(J)
           END DO
        ELSE IF(SUM.NE.ZERO) THEN
           II=I
        END IF
        B(I)=SUM
      END DO

      DO I=N,1,-1
        SUM=B(I)
        DO J=I+1,N
          SUM=SUM-A(I,J)*B(J)
        END DO
        B(I)=SUM/A(I,I)
      END DO

      RETURN

      END

!                                                                      C
!     SUBROUTINE KCOPYMAT                                              C
!                                                                      C
      SUBROUTINE KCOPYMAT(A,B,N)

      IMPLICIT REAL*8(A-H,O-Z)

      DIMENSION A(N,N),B(N,N)

      CALL KCLEAR(B,N,N)

      DO K1=1,N
        DO K2=1,N
          B(K1,K2)=A(K1,K2)
        END DO
      END DO

      RETURN

      END
