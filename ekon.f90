          
      !!!!!! ѕодпрограмма расчета стоимости электроэнергии при заданной цене тепла!!!!!!!!!!!!!!!!!!!
      ! CEL - стоимости электроэнергии, цент/к¬т.
      ! NSTR - число лет строительства,
      ! NEXPL - число лет эксплуатации,
      ! QKBL - капиталовложени€, $,
      ! QNPOLE - полезна€ мощность станции, к¬т,
      ! H - число часов использовани€ установленной мощености,
      ! BCHAS - часовое потребление топлива, тут./час,
      ! STPL - цена газа руб.,
      ! ALYP - дол€ условно-посточнных издержек,
      ! ALAM- дол€ амортизационных издержек,
      ! DIRR- внутренн€€ номра возврата.
      ! ELGODP-годовой расход газа, кг
       ! TEPLGOD-годова€ электрическа€ мощность компрессорных станций
       !STEPL-цена э/э руб к¬т
       !BTOPGOD-годовой расход топлива, ту.т. год
      !!!!  ѕример:
      !!!!    NSTR=4.
      !!!!   NEXPL=30.
      !!!!   ALYP=0.04
      !!!!   ALAM=0.03
      !!!!   DIRR=0.15
      !!!!   STTOP=100. ! долл/ту.т.

      
      SUBROUTINE CELTEP1(CEL,NSTR,NEXPL,QKBL,ELGODP,QKBLST,ALYPST,BTOPGOD,STPL,ALYPLIN,ALAM,DIRR,PCH)
      REAL :: CEL,QKBL,STPL,ALYPLIN,ALAM,DIRR,NSTR,NEXPL,utepl, QKBLLIN
      REAL :: SIGSTR,SIGEX,UT,UYP,EGOD,SIG,ELGODP,QKBLST,ALYPST,BTOPGOD,pch
      integer :: I
      
      STPL=4.2
      ALAM=0.
      ALYPLIN=0.035
      ALYPST=0.06

205   format(/5X,'==== ECONOMY === ')
      206   format(/5X,'QKBL=',F25.8,2X,'SIGSTR=',F16.8,2X,'SIGEX=',F16.8,/5X,'UT=',F25.8,2X,'UYP=',F25.8,/5X, &
       & 'EGOD=',F25.8,2X, 'QKBLST=',F25.8,2X,'SIG=',F16.8,/5X,  'BTOPGOD=',F25.8,2X, 'CEL=',f16.8)
207   format(/5X, 'STPL(gas)=',F16.3,2X,'ALYPST=',F16.3,2X,'DIRR=',f16.3,2x)
      SIGSTR=0.
      DO I=0,int(NSTR-1)
        SIGSTR=SIGSTR+1/((1+DIRR))**I
      END DO
      SIGEX=0.
      DO I=int(NSTR),int(NEXPL)
        SIGEX=SIGEX+1/((1+DIRR))**I
      END DO
      UT=BTOPGOD*STPL
      QKBLLIN=QKBL-QKBLST
      UYP=QKBLLIN*(ALYPLIN+ALAM) + QKBLST*(ALYPST+ALAM)
      EGOD=ELGODP
      
      UTEPL=0.
      
      SIG=SIGSTR/SIGEX/NSTR
      !CEL=((QKBL*SIG+UT+UYP-UTEPL)/EGOD)*100.
      CEL=((QKBL*SIG+UT+UYP)/EGOD)

      
      if (pch==3.) then
        write(*,205)
        write(*,207)STPL,ALYPST,DIRR
        write(*,206)QKBL,SIGSTR,SIGEX,UT,UYP,EGOD,QKBLST,SIG,BTOPGOD,CEL
      end if
      
      if(PCH==4.) then
          open(1,file='tmp_pech.dat',action='write',position='append')
          write(1,205)
          write(1,207)STPL,ALYPST,DIRR
          write(1,206)QKBL,SIGSTR,SIGEX,UT,UYP,EGOD,QKBLST,SIG,BTOPGOD,CEL
          close(1)
      end if
      
      RETURN
      END SUBROUTINE CELTEP1
    

