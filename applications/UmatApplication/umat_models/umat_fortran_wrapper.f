      SUBROUTINE UMAT_WRAPPER(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,STRAN,DSTRAN,
     2 TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,MATERL,NDI,NSHR,NTENS,
     3 NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,CELENT,
     4 DFGRD0,DFGRD1,NOEL,NPT,KSLAY,KSPT,KSTEP,KINC,MATERIALNUMBER)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,STRAN,DSTRAN,
     2 TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,NSHR,
     3 PROPS,COORDS,DROT,PNEWDT,CELENT,
     4 DFGRD0,DFGRD1,NOEL,NPT,KSLAY,KSPT,KSTEP,KINC
C
      INTEGER NTENS,NDI,NSTATV,NPROPS,MATERIALNUMBER,K1,K2
C
      CHARACTER*80 MATERL
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),
     4 DFGRD0(3,3),DFGRD1(3,3)
C
      if(MATERIALNUMBER.EQ.0)THEN
         Call MISES_UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1   RPL,DDSDDT,DRPLDE,DRPLDT,STRAN,DSTRAN,
     2   TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,MATERL,NDI,NSHR,NTENS,
     3   NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,CELENT,
     4   DFGRD0,DFGRD1,NOEL,NPT,KSLAY,KSPT,KSTEP,KINC)
C
      ELSE IF(MATERIALNUMBER.EQ.1)THEN
         Call HYPOPLASTIC_UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1   RPL,DDSDDT,DRPLDE,DRPLDT,STRAN,DSTRAN,
     2   TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,MATERL,NDI,NSHR,NTENS,
     3   NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,CELENT,
     4   DFGRD0,DFGRD1,NOEL,NPT,KSLAY,KSPT,KSTEP,KINC)
C
      ELSE IF(MATERIALNUMBER.EQ.2)THEN
C         Call COMPLIB_NEW(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
C     1   RPL,DDSDDT,DRPLDE,DRPLDT,STRAN,DSTRAN,
C     2   TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,MATERL,NDI,NSHR,NTENS,
C     3   NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,CELENT,
C     4   DFGRD0,DFGRD1,NOEL,NPT,KSLAY,KSPT,KSTEP,KINC)

      ENDIF
C
      RETURN
      END
