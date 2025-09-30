    !Входные
    !sostav gaza: xi - долевой состав газа (sost_dol); mi - молекулярные массы компонентов (m); rost - плотность компонентов
    !tkr0 - критическая температура компонентов;  pkr0 - критическое давление компонентов, МПа
    !ipg - давление газа на входе, кг/см2; itg - температура газа на входе, К; n - число оборотов;
    !igg - расход газа на входе, кг/с; kpompaj - коэффициент помпажирования
    !otg - температура газа на выходе, K; opg - давление газа на выходе, кг/см2
    !chislo_agr - оптимизируемый
    !alpha - вкл выкл определенного типа гпа (4)
    !qrn - высшая объемная теплота сгорания топлива
    !btyr - г.у.т./кВтч удельный расход топлива турбинами 
subroutine gpa_schet6(chislo_agr,ogrrez,q6t3al,qrez1al,q10al,qrez2al,&
    & q12al,qrez3al,q16al,qrez4al,&
    & nnpr,ogr63,ogr10,ogr12,ogr16,ogr,qneb,nnend,kapvl_end,ogg,otg,opg,&
    & igg,ipg,itg,n,kpompaj,btyr,perem_sost,post_sost,alpha,xindr,sost_dol,al6t3,al10,al12,al16,alrez1,alrez2,alrez3,alrez4,pch,hom)
    use raschet_gpa_new
    interface
        subroutine gas_propetries(qtek,qrn,rosm,rsm,z,sost_dol,pin,itg,igg)
            real,intent(in) :: igg,itg,pin
            real, intent(out) ::qtek,qrn,rosm,rsm,z
            real, dimension(7),intent(in) :: sost_dol
        end subroutine gas_propetries
    end interface
        real, dimension (7), intent(in) :: sost_dol
        real, intent(in) :: kpompaj,ipg,igg,btyr,post_sost,pch,hom,xindr
        real :: itg
        
        real, dimension (5), intent(out) :: ogr !(6;6,3;10;12;16)
        real, intent(out) :: opg,otg,qneb,ogg,kapvl_end,nnend
                        
        real, dimension(4), intent(in) :: alpha,n,perem_sost !(6;6,3;10;12;16)
        real, dimension(4),parameter :: zgpa=(/0.933,0.89,0.923,0.9/) !номинальные параметры ГПА 6;6,3;10;12;16МВт
        real, dimension(4),parameter :: rgpa=(/507.,441.4,507.9,507./) !номинальные параметры ГПА 6;6,3;10;12;16МВт
        real, dimension(4),parameter :: tgpa=(/288.,293.,288.,288./) !номинальные параметры ГПА 6;6,3;10;12;16МВт
        real, dimension(4),parameter :: qnom = (/91.,61.,400.,540./) !номинальные параметры ГПА 6;6,3;10;12;16МВт
        real, dimension(4),parameter :: nnom=(/8200.,9000.,6500.,5300./) !номинальные параметры ГПА 6;6,3;10;12;16МВт
        real, dimension(4),parameter :: NGPAnom=(/6300.,10000.,12000.,16000./) !номинальные параметры ГПА 6;6,3;10;12;16МВт
        real, parameter :: kn = 0.95, kobl = 1., ky = 0.985, kt = 3.
        real, parameter :: patm = 0.0968, tcr = 271.75 !атмосферное давление и средняя температура местности (Жигалово)
        real, dimension(4) :: freq,ni,nm,nn,t2,p2,racxodpr,gtyr,Nrasp,Nustan !(6;6,3;10;12;16)
        real, dimension(4), intent(out) :: nnpr
        real, dimension(5), intent(out) :: ogr10,ogr12,ogr16,ogr63
        real, dimension(4) :: sgatie,kpd,racxod_for_one
        real :: pin,qrn,qtek,rosm,rsm,z,gtyr_end,sgatie_end
        !!!=== Новый метод === 
        real, dimension(12) :: al6t3,q6t3al
        real :: gkoef,alrez4,qrez4al
        real, dimension(4) :: nrez,chislo_agr,al16,ogrrez,q16al
        real, dimension(3) :: alrez1,qrez1al
        real, dimension(7) :: al10,q10al
        real, dimension(2) :: alrez2,alrez3,qrez2al,qrez3al
        real, dimension(6) :: al12,q12al
        integer :: i
        !!!===
        !real :: t2end,sgatie_end,p2end,kpd_end
        !!!!!
        !!!=== Новый метод ===
        gkoef=1.
        
        chislo_agr(1)=2. ! не может быть меньше 2 (1 - в работе, 1 - резерв)
        do i=1,12
            chislo_agr(1)=chislo_agr(1) + al6t3(i)
            q6t3al(i)=al6t3(i)*(1.-al6t3(i))*gkoef*alpha(1)
        end do
        nrez=1. ! минимальное число резервных агрегатов любого типа
        do i=1,3
            nrez(1)=nrez(1)+alrez1(i)
            qrez1al(i)=alrez1(i)*(1.-alrez1(i))*gkoef*alpha(1)
        end do
        ogrrez(1)=chislo_agr(1)-1.-nrez(1)
        
        chislo_agr(2)=2.
        do i=1,7
            chislo_agr(2)=chislo_agr(2)+al10(i)
            q10al(i)=al10(i)*(1.-al10(i))*gkoef*alpha(2)
        end do
        nrez(2)=1.
        do i=1,2
            nrez(2)=nrez(2)+alrez2(i)
            qrez2al(i)=alrez2(i)*(1.-alrez2(i))*gkoef*alpha(2)
        end do
        ogrrez(2)=chislo_agr(2)-1.-nrez(2)
        
        chislo_agr(3)=2.
        do i=1,6
            chislo_agr(3)=chislo_agr(3)+al12(i)
            q12al(i)=al12(i)*(1.-al12(i))*gkoef*alpha(3)
        end do
        nrez(3)=1
        do i=1,2
            nrez(3)=nrez(3)+alrez3(i)
            qrez3al(i)=alrez3(i)*(1.-alrez3(i))*gkoef*alpha(3)
        end do
        ogrrez(3)=chislo_agr(3)-1.-nrez(3)
        
        chislo_agr(4)=2.
        do i=1,4
            chislo_agr(4)=chislo_agr(4)+al16(i)
            q16al(i)=al16(i)*(1.-al16(i))*gkoef*alpha(4)
        end do
        nrez(4)=1.
        nrez(4)=nrez(4)+alrez4
        qrez4al=alrez4*(1.-alrez4)*gkoef*alpha(4)
        ogrrez(4)=chislo_agr(4)-1.-nrez(4)
        !!!===
        itg=300.
        !!!!!
        pin=ipg*0.0980665 !давление на входе В МПа
        call gas_propetries(qtek,qrn,rosm,rsm,z,sost_dol,pin,itg,igg)!Находим расход,плотность газа, R, коэфф. сжатия,теплоту сгорания
        freq=n/nnom !отношение кол-ва оборотов к номинальному
        nnpr=freq*SQRT(zgpa*rgpa*tgpa/(rsm*itg*z)) !приведенное отношение оборотов для каждого типа
        racxodpr=qtek*SQRT(zgpa*rgpa*tgpa/(rsm*itg*z)) !приведенный расход газа для каждого типа        
        !racxodpr=qtek*nnom/n !приведенный расход газа для каждого типа - старый вариант, как в литературе 
        !!! старый метод
        !racxod_for_one=racxodpr/(chislo_agr-1.) !расход на 1 ГПА
        !!! ===
        !!! === новый метод ===
        racxod_for_one=racxodpr/(chislo_agr-nrez)
        !!! ===
        if (pch==3.) then
            10 format(5x,10('==='),' ENTER PARAMETERS GPA#',i3,1x,15('==='),//5x,'IGG:',f16.6,2x,'IPG:',f16.6,2x,'PIN:',f16.6,2x,'ITG:',f16.5,&
               & /5x,'Z:',f16.6,4x,'RSM:',f16.6,2x,'ROSM:',f16.6,/5x,'QTEK:',f16.6)
            15 format(/5x,'ALPHA:',4(f16.6),/5x,'N:',4x,4(f16.6),/5x,'FREQ:',1x,4(f16.6),/5x,'NNPR:',1x,4(f16.6))
            25 FORMAT(/5x,'RACXODPR:    ',4(f16.6),/5x,'RACXOD_FOR_1:',4(f16.6))            
            
            
            20 Format(/5x,'CHISLO_AGR:',5x,4(f16.6),/5x,'NREZ:',11x,4(f16.6),/5x,'CHISLO_AGR-NREZ:',4(f16.6))
            
                        
            write(*,10)int(hom),igg,ipg,pin,itg,z,rsm,rosm,qtek
            write(*,15)alpha,n,freq,nnpr
            write(*,20)chislo_agr,nrez,chislo_agr-nrez
            write(*,25)racxodpr,racxod_for_one

        end if

        
        call gpa_NC63(ogr63,kpd(1),sgatie(1),racxod_for_one(1),nnpr(1),kpompaj,pin) !6,3Мвт
        call gpa_108_81_1(ogr10,kpd(2),sgatie(2),racxod_for_one(2),nnpr(2),kpompaj,pin) !10МВт
        call gpa_295_24_1(ogr12,kpd(3),sgatie(3),racxod_for_one(3),nnpr(3),kpompaj,pin) !12МВт
        call gpa_395_21_1(ogr16,kpd(4),sgatie(4),racxod_for_one(4),nnpr(4),kpompaj,pin) !16Мвт
        
        if (pch==3.) then
            90 FORMAT(/5x,'SGATIE:',4(f16.6),/5x,'KPD:   ',4(f16.6))
               write(*,90)sgatie,kpd
        end if
        

        t2=itg*sgatie**0.245 !температура на выходе, К
        p2=pin*sgatie !давление на выходе, МПа

        ni=55.6*pin*racxod_for_one*(sgatie**0.3-1.)/kpd !внутренняя мощность ЦБН 6,10,12,16 МВт
        nm=1-100./NGPAnom !мех. кпд ЦБН 6,10,12,16 МВт

        Nn=ni/(nm*0.95) ! Мощность на муфте ГТУ-ЦБН, 6,10,12,16 кВт одного агрегата
        
        Nrasp=NGPAnom*kn*kobl*ky*(1-kt*(tcr-tgpa)/tcr)*patm/0.1013 !располагаемая мощность ГПА
        
        if (pch==3.) then
            60  FORMAT(/5x,'Ni:   ',4(f16.6),/5x,'Nm:   ',4(f16.6),/5x,'Nn:   ',4(f16.6),/5x,'Nrasp:',4(f16.6))
                write(*,60)ni,nm,Nn,Nrasp
        end if

        
        gtyr=(btyr*7./11.06/3600.)*Nn*(chislo_agr-1.) !расход топлива турбиной 6,10,12,16 МВт
        
        ogr63(5)=Nrasp(1)-(Nn(1)+100.)
        ogr10(5)=Nrasp(2)-(Nn(2)+100.)
        ogr12(5)=Nrasp(3)-(Nn(3)+100.)
        ogr16(5)=Nrasp(4)-(Nn(4)+100.)
        
        if (pch==3.) then
            80  FORMAT(/5x,'OGR63:',5(f16.6),/5x,'OGR10:',5(f16.6),/5x,'OGR12:',5(f16.6),&
               & /5x,'OGR16:',5(f16.6))
                write(*,80)ogr63,ogr10,ogr12,ogr16
        end if
        
		sgatie_end=sum((sgatie-1.)*alpha)+1. !сжатие на выходе

		NNend=sum((Nn+100.)*alpha*chislo_agr) !Установленная мощность ГПА, кВт
		gtyr_end=sum(gtyr*alpha)
		ogg=igg-gtyr_end !расход газа на выходе
        
        opg=(ipg*sgatie_end) !давление газа на выходе ГПА, кг/см2
        otg=sum((t2-itg)*alpha)+itg !температура на выходе ГПА, К
        
        ogr(1)=ogr63(1)*alpha(1)+ogr10(1)*alpha(2)+ogr12(1)*alpha(3)+ogr16(1)*alpha(4)
        ogr(2)=ogr63(2)*alpha(1)+ogr10(2)*alpha(2)+ogr12(2)*alpha(3)+ogr16(2)*alpha(4)
        ogr(3)=ogr63(3)*alpha(1)+ogr10(3)*alpha(2)+ogr12(3)*alpha(3)+ogr16(3)*alpha(4)
        ogr(4)=ogr63(4)*alpha(1)+ogr10(4)*alpha(2)+ogr12(4)*alpha(3)+ogr16(4)*alpha(4)
        ogr(5)=ogr63(5)*alpha(1)+ogr10(5)*alpha(2)+ogr12(5)*alpha(3)+ogr16(5)*alpha(4)
        
        Nustan=NGPAnom*chislo_agr
        kapvl_end=sum(alpha*(post_sost+perem_sost*Nustan))
        qneb=xindr-sum(alpha)
        if (pch==3.) then
               40 FORMAT(/5x,10('==='),' OUTPUT PARAMETERS GPA#',i3,1x,'=',21('=='),//5x,'OGG:',f16.6,2x,'OPG:',f16.6,2x,'OTG:',f16.6,&
               & /5x,'NPOTR:',f16.6,2x,'NUSTAN:',f16.6,/5x,'GTYR:',f16.6,2x,'SGATIE:',f16.6,//5x,'OGR:',5(f16.6),//5x,'KAPVL:',f30.6)
               110 FORMAT(/5x,'XINDR:',f16.6,2x,'QNEB:',f16.6,//5x,100('='))
               150 format(/5x,50('--'),/5x,'AL6t3:',2x,6(f16.6),/13x,6(f16.6),/5x,'Q6t3AL:',1x,6(f16.6),/13x,6(f16.6),/5x,'ALREZ1:',1x,3(f16.6),/5x,'QREZ1AL:',3(f16.6),&
               & //5x,'AL10:',3x,6(f16.6),/13x,f16.6,/5x,'Q10AL:',2x,6(f16.6),/13x,f16.6,/5x,'ALREZ2:',1x,2(f16.6),/5x,'QREZ2AL:',2(f16.6),&
               & //5x,'AL12:',3x,6(f16.6),/5x,'Q12AL:',2x,6(f16.6),/5x,'ALREZ3:',1x,2(f16.6),/5x,'QREZ3AL:',2(f16.6),&
               & //5x,'AL16:',3x,4(f16.6),/5x,'Q16AL:',2x,4(f16.6),/5x,'ALREZ4:',1x,f16.6,/5x,'QREZ4AL:',f16.6,&
               & /5x,'OGRREZ:',1x,4(f16.6),/5x,100('-'),/5x,100('='))
            write(*,40)int(hom),ogg,opg,otg,NNend,sum(Nustan*alpha),gtyr_end,sgatie_end,ogr,kapvl_end
            write(*,110)xindr,qneb
            write(*,150)al6t3(1:6),al6t3(7:12),q6t3al(1:6),q6t3al(1:6),alrez1,qrez1al,al10(1:6),al10(7),q10al(1:6),q10al(7),alrez2,qrez2al,al12,q12al,alrez3,qrez3al,al16,q16al,alrez4,qrez4al,ogrrez
        end if
        

end subroutine gpa_schet6