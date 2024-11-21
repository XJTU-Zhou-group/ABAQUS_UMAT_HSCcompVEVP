C********************************TOP***************************
C     UMAT FOR ABAQUS/STANDARD
C     VISCOELASTO-VISCOPLASTIC MODEL FOR CFRP WOVEN-COMPOSITES
C     2D PLANE STRESS
C     THIS CODE WAS DEVELOPED BY XIAOWEI YUE
C*******************************************************************
      
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
      
      INCLUDE 'ABA_PARAM.INC'
      
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)

      PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0,THREE=3.D0,
     1 SIX=6.D0, NINE=9.D0, TOLER=1.D-8)
      
      INTEGER NProny,counter
      REAL*8 rhoi(9)
      REAL*8 C11inf,C12inf,C22inf,C66inf
      REAL*8 C11i(9),C12i(9),C22i(9),C66i(9),C66RM,C12RM,C22RM
      REAL*8 qold(NTENS,9),q(NTENS,9)
      REAL*8 f11,f22,f66
      REAL*8 depsilon(NTENS),epsilonE(NTENS),sigmaE(NTENS)
      REAL*8 STRESSOLD(NTENS),DPSTRAN(NTENS),DSTRESS(NTENS)
      REAL*8 DESTRAN(NTENS),PFDIR(NTENS)
      REAL*8 HISTR,A66,R1,R2,P,YF,XYF,XRES,XPHI,XPHIDP,XPHIR
      
      DO K=1,NTENS
          STRESSOLD(K)=STRESS(K)
      END DO  
      
      
      NProny=9
      DO I=1,NProny
          rhoi(I)=PROPS(I)
          C11i(I)=PROPS(10+I)
          C12i(I)=PROPS(20+I)
          C22i(I)=PROPS(30+I)
          C66i(I)=PROPS(40+I)
      END DO
      
      C11inf=PROPS(10)
      C12inf=PROPS(20)
      C22inf=PROPS(30)
      C66inf=PROPS(40)
      
      A66=PROPS(50)
      SIGY0=PROPS(51)
      XB1=PROPS(52)
      XB2=PROPS(53)
      XQ1=PROPS(54)
      XQ2=PROPS(55)
      XM=PROPS(56)
      XK=PROPS(57)
      
      counter=1
      
      DO J=1,NTENS
          DO I=1,NProny
              qold(J,I)=STATEV(counter)
              counter=counter+1
          END DO
      END DO
      
      
      DO I=1,NTENS
          DO J=1,NTENS
              DDSDDE(I,J)=0.0
          END DO
      END DO
      
      DO I=1,NProny
          DDSDDE(1,1)=DDSDDE(1,1)+
     1    C11i(I)*(1-EXP(-DTIME/rhoi(I)))/(DTIME/rhoi(I))
          DDSDDE(1,2)=DDSDDE(1,2)+
     1    C12i(I)*(1-EXP(-DTIME/rhoi(I)))/(DTIME/rhoi(I))
          DDSDDE(2,2)=DDSDDE(2,2)+
     1    C22i(I)*(1-EXP(-DTIME/rhoi(I)))/(DTIME/rhoi(I))
          DDSDDE(3,3)=DDSDDE(3,3)+
     1    C66i(I)*(1-EXP(-DTIME/rhoi(I)))/(DTIME/rhoi(I))
      END DO
      

      DDSDDE(1,1)=DDSDDE(1,1)+C11inf
      DDSDDE(1,2)=DDSDDE(1,2)+C12inf
      DDSDDE(2,2)=DDSDDE(2,2)+C22inf
      DDSDDE(3,3)=DDSDDE(3,3)+C66inf

      
      DDSDDE(2,1)=DDSDDE(1,2)
      
      C66RM=DDSDDE(3,3)
      C12RM=DDSDDE(1,2)
      C22RM=DDSDDE(2,2)
      
      depsilon(1)=DSTRAN(1)
      depsilon(2)=DSTRAN(2)
      depsilon(3)=DSTRAN(3)
      epsilonE(1)=STRAN(1)+DSTRAN(1)
      epsilonE(2)=STRAN(2)+DSTRAN(2)
      epsilonE(3)=STRAN(3)+DSTRAN(3)
      
      

      DO J=1,NTENS
          DO I=1,NProny
              q(J,I)=EXP(-DTIME/rhoi(I))*qold(J,I)+
     1        DSTRAN(J)*(1-EXP(-DTIME/rhoi(I)))/(DTIME/rhoi(I))
          END DO
      END DO

      
      f11=0
      f22=0
      f66=0
      
      DO I=1,NProny
          f11=f11+(EXP(-DTIME/rhoi(I))-ONE)*
     1    (qold(1,I)*C11i(I)+qold(2,I)*C12i(I))
          f22=f22+(EXP(-DTIME/rhoi(I))-ONE)*
     1    (qold(1,I)*C12i(I)+qold(2,I)*C22i(I))
          f66=f66+(EXP(-DTIME/rhoi(I))-ONE)*qold(3,I)*C66i(I)
      END DO
      
      
      DO K1 = 1,NTENS
         DO K2 = 1,NTENS
              STRESS(K1)=STRESS(K1)+DDSDDE(K1,K2)*DSTRAN(K2)
          END DO
      END DO
      
      STRESS(1)=STRESS(1)+f11
      STRESS(2)=STRESS(2)+f22
      STRESS(3)=STRESS(3)+f66
      

      R1=STATEV(28)
      R2=STATEV(29)
      P=STATEV(30)     
      
      
      HISTR=0.0D0
      HISTR=DSQRT((THREE/TWO)*(TWO*
     1 STRESS(3)*STRESS(3)))
      
      !IF(HISTR.LT.1.D-15)THEN
      !HISTR=1.D-15
      !END IF      
      
      PFDIR(1)=ZERO
      PFDIR(2)=ZERO
      PFDIR(3)=THREE*STRESS(3)/HISTR
      
      
      STATEV(31)=HISTR
      
      YF=HISTR - (R1+R2) - SIGY0
      
      IF(YF.GT.0.)THEN
          R1T=R1
          R2T=R2
          P0=P
          XDP=0.
          DO KNEW=1,10
              XYF=HISTR-THREE*C66RM*XDP-(R1+R2)-SIGY0
              XPHI=(XYF/XK)**(ONE/XM)
              XPHIR=(-ONE/(XM*XK))*((XYF/XK)**((ONE/XM)-ONE))
              XPHIDP=XPHIR*THREE*C66RM  
              XRES=XPHI-(XDP/DTIME)
              DEQPL=XRES/((ONE/DTIME)-XPHIDP-XPHIR*(
     1        XB1*(XQ1-R1)+XB2*(XQ2-R2)))
              XDP=XDP+DEQPL
             ! STATEV(KNEW+7)=XRES
              R1=R1T+XB1*(XQ1-R1T)*XDP
              R2=R2T+XB2*(XQ2-R2T)*XDP
              IF(ABS(XRES).LT.TOLER)GOTO 10
          END DO
10        CONTINUE
      P=P0+XDP    
      DPSTRAN(1)=ZERO
      DPSTRAN(2)=ZERO
      DPSTRAN(3)=THREE*XDP*STRESS(3)/HISTR
      
      STRESS(1)=STRESS(1)
      STRESS(2)=STRESS(2)
      STRESS(3)=STRESS(3)-C66RM*PFDIR(3)*XDP
      
      DESTRAN(1)=DSTRAN(1)-DPSTRAN(1)
      DESTRAN(2)=DSTRAN(2)-DPSTRAN(2)
      DESTRAN(3)=DSTRAN(3)-DPSTRAN(3)      
      
      STATEV(32)=DPSTRAN(1)
      STATEV(33)=DPSTRAN(2)
      STATEV(34)=DPSTRAN(3)
      STATEV(35)=DESTRAN(1)
      STATEV(36)=DESTRAN(2)
      STATEV(37)=DESTRAN(3)
      STATEV(38)=DSTRAN(1)
      STATEV(39)=DSTRAN(2)
      STATEV(40)=DSTRAN(3)
      STATEV(41)=XDP
      
      DO J=1,NTENS
          DO I=1,NProny
              q(J,I)=EXP(-DTIME/rhoi(I))*qold(J,I)+
     1        DESTRAN(J)*(1-EXP(-DTIME/rhoi(I)))/(DTIME/rhoi(I))
          END DO
      END DO

   
      
      !DO K1 = 1,NTENS
      !   DO K2 = 1,NTENS
      !        STRESS(K1)=STRESSOLD(K1)+DDSDDE(K1,K2)*DESTRAN(K2)
      !    END DO
      !END DO
      !
      !STRESS(1)=STRESS(1)+f11
      !STRESS(2)=STRESS(2)+f22
      !STRESS(3)=STRESS(3)+f66
      
      STATEV(28)=R1
      STATEV(29)=R2
      STATEV(30)=P
      
      
      
      END IF      
      

      
      counter=1
      DO J=1,NTENS
          DO I=1,NProny
              STATEV(counter)=q(J,I)
              counter=counter+1
          END DO
      END DO
      


      
      SSE=0.0
      SCD=0.0
      SPD=0.0
      
      
      
      RETURN
      END
