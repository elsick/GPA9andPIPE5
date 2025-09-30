module raschet_gpa_new
contains
    subroutine gpa_108_81_1(q,kpd,sgatie,racxod,freq,kpompaj,pin)
    !10МВт
    !implicit none
    real, parameter :: pmin = 4.086431055, pmax = 12.2583125
    integer,parameter :: nper_sgatie=2, nkpol_sgatie=28, nkpol_sgatie_max=7, nkpol_sgatie_min=7, nkpol_zapiranie=7, nkpol_pompaj=7
    integer,parameter :: nper_kpd=2, nkpol_kpd=21, nkpol_sgatie_min_zapiranie=7, nkpol_sgatie_max_pompaj=7
    real,intent(out) :: sgatie,kpd,q(4)
    real :: zapiranie,pompaj,polinom_one_per !sgatie_max,sgatie_min,
    real,intent(in) :: racxod,freq,kpompaj,pin
    !степень сжатия от расхода и частоты
    real, dimension (nkpol_sgatie),parameter :: kpol_sgatie = (/ 245.304643650698, -1661.86422610218, 4526.50688155926,&
                                                     & -5550.46419562822, 2039.91101199152, -2320.35726028259,&
                                                     & 887.422907262760, 0.433422362796477, 2.48843871148831,&
                                                     & -51.8422754094591, 160.037519416181, 89.3714832948233,&
                                                     & -40.4608094133869, -3.640845975565193E-002, 0.711369509306107,&
                                                     & -2.41817034716001, -4.88012699760684, 0.840813232835174,&
                                                     & -2.800967206702394E-003, 1.307457337206214E-002, 9.374831354654634E-002,&
                                                     & 6.946188998735063E-003, -1.623765324042589E-005, -7.668532629829249E-004,&
                                                     & -4.100054203612867E-004, 2.278307612384356E-006, 4.537668799652108E-006,&
                                                     & -1.617071114481301E-008/)
    integer,  dimension(nkpol_sgatie,nper_sgatie), parameter  :: mst_sgatie = (/ 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2,&
                                                                               & 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 5, 5, 6,&
                                                                               & 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 0,&
                                                                               & 1, 2, 3, 4, 0, 1, 2, 3, 0, 1, 2, 0, 1, 0/)
    !максимальная степень сжатия от расхода
    real, dimension(nkpol_sgatie_max), parameter :: kpol_sgatie_max = (/ -356.743769067928, 30.0785658963337, -1.04190783340311,&
                          & 1.916896193632236E-002, -1.973124016873648E-004, 1.076543914900617E-006, -2.433879034218042E-009/)
    integer,  dimension(nkpol_sgatie_max), parameter  :: mst_sgatie_max = (/ 0, 1, 2, 3, 4, 5, 6/)
    !минимальная степень сжатия от расхода
    real, dimension(nkpol_sgatie_min), parameter :: kpol_sgatie_min = (/-606.599104641618,81.3099560444396,-4.50714624786421,&
        & 0.132641554691307,-2.185557125022853E-003,1.911681831579878E-005,-6.935924011825390E-008/)
    integer,  dimension(nkpol_sgatie_min), parameter  :: mst_sgatie_min =(/0, 1, 2, 3, 4, 5, 6/)
    !линия запирания: расход от частоты
    real, dimension(nkpol_zapiranie), parameter :: kpol_zapiranie = (/-2691.96924273906, 18764.2685636099, -53781.6495463471,&
                          & 82122.4272708306, -70122.4224829242, 31737.3310840753, -5947.49785670401/)
    integer,  dimension(nkpol_zapiranie), parameter  :: mst_zapiranie = (/ 0, 1, 2, 3, 4, 5, 6/)
    !линия помпажа: расход от частоты
    real, dimension(nkpol_pompaj), parameter :: kpol_pompaj = (/ 2449.16774409988, -16736.0338260356, 47731.1029694251,&
                          & -72044.9801491405, 60881.1367304423, -27307.4868041446, 5078.91745038244/)
    integer,  dimension(nkpol_pompaj), parameter  :: mst_pompaj = (/ 0, 1, 2, 3, 4, 5, 6/)
    !к.п.д. от расхода и частоты
    real, dimension(nkpol_kpd), parameter :: kpol_kpd = (/ -1.36027466686898, 7.43992882371781, -30.3144529749605,&
                          & 159.716864373941, -213.803310836350, 2.46069192889359,&
                          & 7.225395444121546E-002, 0.400906956928256, -6.06718413984175,&
                          & 8.65102541284487, 2.82703660743059, -5.347053468306997E-003,&
                          & 8.570359218640838E-002, -0.113366270708348, -0.146707434871160,&
                          & -3.583388741965287E-004, 3.426416118499013E-004, 2.740518870738017E-003,&
                          & 1.178657147072004E-006, -2.165724090152915E-005, 6.097294571613177E-008/)
   integer,  dimension(nkpol_kpd,nper_kpd), parameter  :: mst_kpd = (/ 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 4, 5, 0, 1, 2,&
                          & 3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 3, 0, 1, 2, 0, 1, 0/)
   !линия запирания, степень сжатия от расхода
    real, dimension(nkpol_sgatie_min_zapiranie), parameter :: kpol_sgatie_min_zapiranie = (/-1198.73014636578,101.360310601609,-3.55148110195914,&
        & 6.606957252332668E-002,-6.882335151815197E-004,3.806789236330611E-006,-8.735021414704294E-009/)
    integer, dimension(nkpol_sgatie_min_zapiranie) :: mst_sgatie_min_zapiranie = (/0,1,2,3,4,5,6/)
    !линия помпажа, степень сжатия от расхода
    real,dimension(nkpol_sgatie_max_pompaj), parameter :: kpol_sgatie_max_pompaj = (/1073.91830224567,-142.603947661004,7.85623356668443,&
        & -0.229604284676414,3.755645980713423E-003,-3.259692357672826E-005,1.173118685766874E-007/)
    integer, dimension(nkpol_sgatie_max_pompaj) :: mst_sgatie_max_pompaj = (/0,1,2,3,4,5,6/)
    !!!РАСЧЕТ!!!
    !от частоты находим расход по линии запирания
    zapiranie=polinom_one_per(kpol_zapiranie,mst_zapiranie,nkpol_zapiranie,freq)
    !от частоты находим расход по линии помпажа
    pompaj=polinom_one_per(kpol_pompaj,mst_pompaj,nkpol_pompaj,freq)
    q(1)=racxod-pompaj*kpompaj
    q(2)=zapiranie-racxod
    q(3)=pin-pmin
    if ((q(1)>=0.).and.(q(2)>=0.).and.(q(3) >=0.)) then
        !от частоты и расхода находим степень сжатия
        call Polinom (sgatie,kpol_sgatie,mst_sgatie,nkpol_sgatie,nper_sgatie,racxod,freq)
        !от частоты и расхода находим к.п.д.
        call Polinom(kpd,kpol_kpd,mst_kpd,nkpol_kpd,nper_kpd,racxod,freq)
    else
        sgatie=1.
        kpd=1.
    end if
    !!от расхода находим максимальное сжатие
    !if ((56.903 < racxod).and.(racxod <=88.576)) then
    !    sgatie_max=polinom_one_per(kpol_sgatie_max,mst_sgatie_max,nkpol_sgatie_max,racxod)
    !else if ((36.255 < racxod).and.(racxod <= 56.903)) then
    !    sgatie_max=polinom_one_per(kpol_sgatie_max_pompaj,mst_sgatie_max_pompaj,nkpol_sgatie_max_pompaj,racxod)
    !else 
    !    sgatie_max=3.957
    !end if
    !!от расхода находим минимальное сжатие
    !if ((56.134 < racxod).and.(racxod <= 88.576)) then
    !    sgatie_min=polinom_one_per(kpol_sgatie_min_zapiranie,mst_sgatie_min_zapiranie,nkpol_sgatie_min_zapiranie,racxod)
    !else if ((36.255 < racxod).and.(racxod <= 56.134)) then
    !    sgatie_min=polinom_one_per(kpol_sgatie_min,mst_sgatie_min,nkpol_sgatie_min,racxod)
    !else 
    !    sgatie_min=1.575
    !end if


    
    !Ограничения
    q(4)=pmax-pin*sgatie

    end subroutine gpa_108_81_1
    
 subroutine gpa_295_24_1(q,kpd,sgatie,racxod,freq,kpompaj,pin)
    !12МВт
    real, parameter :: pmin = 3.81478685, pmax = 5.491724
    integer,parameter :: nper_sgatie=2, nkpol_sgatie=28, nkpol_sgatie_max=7, nkpol_sgatie_min=7, nkpol_zapiranie=7, nkpol_pompaj=7
    integer,parameter :: nper_kpd=2, nkpol_kpd=28, nkpol_sgatie_min_zapiranie=7, nkpol_sgatie_max_pompaj=7
    real,intent(out) :: sgatie,kpd,q(4)
    real :: zapiranie,pompaj,polinom_one_per !sgatie_max,sgatie_min,
    real,intent(in) :: racxod,freq,kpompaj,pin
    !степень сжатия от расхода и частоты
    real, dimension (nkpol_sgatie),parameter :: kpol_sgatie = (/ -17.8190037670543, 142.523096684284, -448.258258634105, &
    & 746.526236978471, -683.627923124701, 317.859673372094, -51.3469741427598, -3.008384868895519E-002, 0.204918887032357, &
    & -0.507471871510986, 0.510723629195493, -0.103953810207499, -0.145026724547944, -5.175810702361068E-005, 1.508221563736929E-004, &
    & 2.322672655322826E-004, -1.051433495159146E-003, 1.184007606751029E-003, 6.012409380865421E-008, -8.758028464612733E-007, &
    & 2.526141887785887E-006, -3.281424481618426E-006, 4.369012821395046E-010, -2.247520154205496E-009, 4.756798816445988E-009, &
    & 7.691768377902528E-013, -3.715928113409995E-012, 1.231725891236637E-015/)
    integer,  dimension(nkpol_sgatie,nper_sgatie), parameter  :: mst_sgatie = (/ 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2,&
                                                                               & 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 5, 5, 6,&
                                                                               & 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 0,&
                                                                               & 1, 2, 3, 4, 0, 1, 2, 3, 0, 1, 2, 0, 1, 0/)
    !максимальная степень сжатия от расхода
    real, dimension(nkpol_sgatie_max), parameter :: kpol_sgatie_max = (/ 0.969809413933750, 8.537679351854503E-003, -5.536664625133838E-005, &
    & 2.099917494908951E-007, -4.264890175349893E-010, 4.074326950702102E-013, -1.427048621138134E-016/)
    integer,  dimension(nkpol_sgatie_max), parameter  :: mst_sgatie_max = (/ 0, 1, 2, 3, 4, 5, 6/)
    !минимальная степень сжатия от расхода
    real, dimension(nkpol_sgatie_min), parameter :: kpol_sgatie_min = (/-1.81474718831348, 6.769653234249770E-002, -6.279402990235789E-004, &
    & 3.101324007484858E-006, -8.525182483100763E-009, 1.225934347794058E-011, -7.200536647176040E-015/)
    integer,  dimension(nkpol_sgatie_min), parameter  :: mst_sgatie_min =(/0, 1, 2, 3, 4, 5, 6/)
    !линия запирания: расход от частоты
    real, dimension(nkpol_zapiranie), parameter :: kpol_zapiranie = (/46300.7191958765, -311288.508392047, 870226.814052511, &
    & -1288155.11905359, 1066745.22419829, -468617.682870566, 85325.1772657454/)
    integer,  dimension(nkpol_zapiranie), parameter  :: mst_zapiranie = (/ 0, 1, 2, 3, 4, 5, 6/)
    !линия помпажа: расход от частоты
    real, dimension(nkpol_pompaj), parameter :: kpol_pompaj = (/ -36273.9080051499, 249477.812756421, -709225.176361105, &
    &1070099.75218802, -902960.247777459, 404060.530753240, -74922.6503838445/)
    integer,  dimension(nkpol_pompaj), parameter  :: mst_pompaj = (/ 0, 1, 2, 3, 4, 5, 6/)
    !к.п.д. от расхода и частоты
    real, dimension(nkpol_kpd), parameter :: kpol_kpd = (/ 55.5675036137725, -357.333197211234, 945.636189011430, &
    & -1325.60751658845, 1056.50306577784, -381.981951085169, -51.8461568145034, -3.516458531304605E-002, 0.309790269150171, &
    & -0.769017130492720, 0.707390129994302, -1.27052514511824, 2.03869268661652, -1.526554651709129E-004, 4.589543463802873E-004, &
    & 9.955774452580216E-005, 4.937256798590817E-003, -1.198856529967827E-002, -1.833128412353176E-008, -1.035006969826144E-006, &
    & -1.321216498991680E-005, 3.799820686469392E-005, 6.290897712739772E-010, 1.821600191166202E-008, -6.587670827574631E-008, &
    & -9.762748593502647E-012, 5.876705485701466E-011, -2.098767359367412E-014/)
    integer,  dimension(nkpol_kpd,nper_kpd), parameter  :: mst_kpd = (/ 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, &
    & 4, 4, 4, 5, 5, 6, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 3, 0, 1, 2, 0, 1, 0/)
    !линия запирания, степень сжатия от расхода
    real, dimension(nkpol_sgatie_min_zapiranie), parameter :: kpol_sgatie_min_zapiranie = (/-32.7206662700630,0.427438917520890,-2.247001768006278E-003,&
    & 6.274584425335627E-006,-9.802813876539310E-009,8.124887668193730E-012,-2.790208138503165E-015/)
    integer,  dimension(nkpol_sgatie_min_zapiranie), parameter  :: mst_sgatie_min_zapiranie = (/0, 1, 2, 3, 4, 5, 6/)
    !линия помпажа, степень сжатия от расхода
    real, dimension(nkpol_sgatie_max_pompaj), parameter :: kpol_sgatie_max_pompaj = (/-36.9699845035810,0.899362829636555,-8.691456674876521E-003,&
    & 4.381666205681122E-005,-1.206065798426946E-007,1.706149334790116E-010,-9.544773640841157E-014/)
    integer,  dimension(nkpol_sgatie_max_pompaj), parameter  :: mst_sgatie_max_pompaj = (/0, 1, 2, 3, 4, 5, 6/)
    !!!РАСЧЕТ!!!
    !от частоты находим расход по линии запирания
    zapiranie=polinom_one_per(kpol_zapiranie,mst_zapiranie,nkpol_zapiranie,freq)
    !от частоты находим расход по линии помпажа
    pompaj=polinom_one_per(kpol_pompaj,mst_pompaj,nkpol_pompaj,freq)
    q(1)=racxod-pompaj*kpompaj
    q(2)=zapiranie-racxod
    q(3)=pin-pmin
    if ((q(1) >= 0.).and.(q(2) >= 0.).and.(q(3) >= 0.)) then
        !от частоты и расхода находим степень сжатия
        call Polinom (sgatie,kpol_sgatie,mst_sgatie,nkpol_sgatie,nper_sgatie,racxod,freq)
        !от частоты и расхода находим к.п.д.
        call Polinom(kpd,kpol_kpd,mst_kpd,nkpol_kpd,nper_kpd,racxod,freq)
    else
        sgatie=1.
        kpd=1.
    end if
    
    q(4)=pmax-pin*sgatie

    
end subroutine gpa_295_24_1
    
subroutine gpa_395_21_1(q,kpd,sgatie,racxod,freq,kpompaj,pin)
    !16 МВт
    real, parameter :: pmin = 3.87362675, pmax = 7.453054
    integer,parameter :: nper_sgatie=2, nkpol_sgatie=28, nkpol_sgatie_max=7, nkpol_sgatie_min=7, nkpol_zapiranie=8, nkpol_pompaj=7
    integer,parameter :: nper_kpd=2, nkpol_kpd=28, nkpol_sgatie_min_zapiranie=7, nkpol_sgatie_max_pompaj=7
    real,intent(out) :: sgatie,kpd,q(4)
    real :: sgatie_max,sgatie_min,zapiranie,pompaj,polinom_one_per
    real,intent(in) :: racxod,freq,kpompaj,pin
    !степень сжатия от расхода и частоты
    real, dimension (nkpol_sgatie),parameter :: kpol_sgatie = (/ 15.3804766271747, -107.398570558461, 346.201888832740, &
    & -612.386055437106, 627.708197925773, -346.241754535994, 77.6842906642766, 3.129990991456906E-002, -0.253249037856440, &
    & 0.798233398737456, -1.26315347275004, 0.927925098831265, -0.234145386584284, 8.619904884785587E-005, -5.362942809816210E-004, &
    & 1.399201397340330E-003, -1.204817820229740E-003, 2.153748571702565E-004, 1.173917293572982E-007, -8.985379261737802E-007, &
    & 5.875478942784972E-007, 3.384592539710540E-007, 3.367451334650449E-010, 3.749374850968549E-010, &
    & -1.024855285768825E-009, -4.998951759677232E-013, 8.443947336976687E-013, -1.485830804003981E-016/)
    integer,  dimension(nkpol_sgatie,nper_sgatie), parameter  :: mst_sgatie = (/ 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2,&
                                                                               & 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 5, 5, 6,&
                                                                               & 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 0,&
                                                                               & 1, 2, 3, 4, 0, 1, 2, 3, 0, 1, 2, 0, 1, 0/)
    !максимальная степень сжатия от расхода
    real, dimension(nkpol_sgatie_max), parameter :: kpol_sgatie_max = (/ 1.99628383405842, -8.051131073965555E-003, 5.739366659520784E-005, &
    & -1.923667180714285E-007, 3.430606406866368E-010, -3.351014668557940E-013, 1.414110617480467E-016/)
    integer,  dimension(nkpol_sgatie_max), parameter  :: mst_sgatie_max = (/ 0, 1, 2, 3, 4, 5, 6/)
    !минимальная степень сжатия от расхода
    real, dimension(nkpol_sgatie_min), parameter :: kpol_sgatie_min = (/-31.6029484647928, 0.734909290390101, -6.781227332015312E-003, &
    & 3.301667835446129E-005, -8.944478597732114E-008, 1.277870095990796E-010, -7.523196002917723E-014/)
    integer,  dimension(nkpol_sgatie_min), parameter  :: mst_sgatie_min =(/0, 1, 2, 3, 4, 5, 6/)
    !линия запирания: расход от частоты
    real, dimension(nkpol_zapiranie), parameter :: kpol_zapiranie = (/-339966.908946367, 2799134.86165256, -9816814.13889742, &
    & 19023652.2176249, -21998580.9611296, 15182163.1767442, -5790926.08568809, 941877.014490333/)
    integer,  dimension(nkpol_zapiranie), parameter  :: mst_zapiranie = (/ 0, 1, 2, 3, 4, 5, 6, 7/)
    !линия помпажа: расход от частоты
    real, dimension(nkpol_pompaj), parameter :: kpol_pompaj = (/ -47902.7094068059, 324380.758975730, -909415.572231341, &
    & 1354685.44035731, -1130061.85267271, 500603.905586597, -92018.1226327659/)
    integer,  dimension(nkpol_pompaj), parameter  :: mst_pompaj = (/ 0, 1, 2, 3, 4, 5, 6/)
    !к.п.д. от расхода и частоты
    real, dimension(nkpol_kpd), parameter :: kpol_kpd = (/ 4.95999156406392, -48.5359085587714, 174.768739489912, &
    & -280.569274524175, 225.015394685675, -123.194708282897, 292.118164237543, 5.381687851006510E-002, -0.208798231713009, &
    & 0.123652082682579, 0.211219653473148, 0.204089655306217, -4.32071636378986, -1.100865366666807E-004, 9.916640772846981E-004, &
    & -1.750901938055132E-003, -1.285703929722736E-003, 2.832240687063669E-002, -5.339604605313858E-007, 7.399543620559355E-007, &
    & 6.000222646084568E-006, -9.787562631503782E-005, 3.868120803533278E-010, -8.301182642136284E-009, &
    & 1.864636529043660E-007, 3.736379839166495E-012, -1.870753568149834E-010, 7.753105729791131E-014/)
   integer,  dimension(nkpol_kpd,nper_kpd), parameter  :: mst_kpd = (/ 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, &
    & 4, 4, 4, 5, 5, 6, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 3, 0, 1, 2, 0, 1, 0/)
    !Линия запирания, степень сжатия от расхода
    real, dimension(nkpol_sgatie_min_zapiranie), parameter :: kpol_sgatie_min_zapiranie = (/-177.539304305353,2.23985697186363,-1.165987949955057E-002,&
    & 3.224352484866068E-005,-4.994606347086464E-008,4.109357199466589E-011,-1.402939894572560E-014/)
    integer,  dimension(nkpol_sgatie_min_zapiranie), parameter  :: mst_sgatie_min_zapiranie =(/0, 1, 2, 3, 4, 5, 6/)
    !Линия помпажа степень сжатия от расхода
    real, dimension(nkpol_sgatie_max_pompaj), parameter :: kpol_sgatie_max_pompaj = (/ 73.1486336996584,-1.81986603842921,1.900545878218060E-002,&
    & -1.050968124148446E-004,3.249045118748617E-007,-5.322746280470412E-010,3.611219762031987E-013/)
    integer,  dimension(nkpol_sgatie_max_pompaj), parameter  :: mst_sgatie_max_pompaj = (/ 0, 1, 2, 3, 4, 5, 6/)
    !!!РАСЧЕТ!!!
    !от частоты находим расход по линии запирания
    zapiranie=polinom_one_per(kpol_zapiranie,mst_zapiranie,nkpol_zapiranie,freq)
    !от частоты находим расход по линии помпажа
    pompaj=polinom_one_per(kpol_pompaj,mst_pompaj,nkpol_pompaj,freq)
    q(1)=racxod-pompaj*kpompaj
    q(2)=zapiranie-racxod
    q(3)=pin-pmin
    if ((q(1) >= 0.).and.(q(2) >= 0.).and.(q(3) >= 0.)) then
        !от частоты и расхода находим степень сжатия
        call Polinom (sgatie,kpol_sgatie,mst_sgatie,nkpol_sgatie,nper_sgatie,racxod,freq)
        !от частоты и расхода находим к.п.д.
        call Polinom(kpd,kpol_kpd,mst_kpd,nkpol_kpd,nper_kpd,racxod,freq)
    else
        sgatie=1.
        kpd=1.
    end if
    !от расхода находим максимальное сжатие
    !if ((299.201 <racxod).and.(racxod <= 593.485)) then
    !    sgatie_max=polinom_one_per(kpol_sgatie_max,mst_sgatie_max,nkpol_sgatie_max,racxod)
    !else if ((190.897 <racxod).and.(racxod <= 299.201)) then
    !    sgatie_max=polinom_one_per(kpol_sgatie_max_pompaj,mst_sgatie_max_pompaj,nkpol_sgatie_max_pompaj,racxod)
    !else
    !    sgatie_max=1.619
    !end if
    !от расхода находим минимальное сжатие
    !if ((190.987 < racxod).and.(racxod <= 380.949)) then
    !    sgatie_min=polinom_one_per(kpol_sgatie_min,mst_sgatie_min,nkpol_sgatie_min,racxod)
    !else if ((380.949 < racxod).and.(racxod <= 593.485)) then
    !    sgatie_min=polinom_one_per(kpol_sgatie_min_zapiranie,mst_sgatie_min_zapiranie,nkpol_sgatie_min_zapiranie,racxod)
    !else
    !    sgatie_min=1.108
    !end if
    !Ограничения
    !q(1)=sgatie_max-sgatie
    !q(2)=sgatie-sgatie_min
    
    q(4)=pmax-pin*sgatie

end subroutine gpa_395_21_1
    
subroutine gpa_2H637(q,kpd,sgatie,racxod,freq,kpompaj,pin)
    !6 МВт 2Н-6-37-1,6
    real, parameter :: pmin = 2.06920315, pmax = 3.6284605
    integer,parameter :: nper_sgatie=2, nkpol_sgatie=28, nkpol_sgatie_max=9, nkpol_sgatie_min=9, nkpol_zapiranie=7, nkpol_pompaj=7
    integer,parameter :: nper_kpd=2, nkpol_kpd=28, nkpol_sgatie_min_zapiranie=7, nkpol_sgatie_max_pompaj = 7
    real,intent(out) :: sgatie,kpd,q(4)
    real :: sgatie_max,sgatie_min,zapiranie,pompaj,polinom_one_per
    real,intent(in) :: racxod,freq,kpompaj,pin
    !степень сжатия от расхода и частоты
    real, dimension (nkpol_sgatie),parameter :: kpol_sgatie = (/184.735769352106,-1353.33440926842,4048.27140600656,&
        & -6231.28812238778,5144.14302242215,-2093.48052387895,&
        & 326.825882118426,0.400363163321675,-1.64076180797726,0.896516137928745,3.94218941326244,-6.82632836328563,&
        & 2.37462125268623,-1.519712709527136E-003,1.476829940418137E-002,-4.073190162257383E-002,5.290885482635060E-002,&
        & -1.290257740345910E-002,-1.486027831697334E-005,6.137652283511266E-005,-1.347555890535293E-004,-9.190338240235045E-006,&
        & -1.873059424419826E-008,2.012541794932363E-007,2.415293248068948E-007,-1.893823631169773E-010,-7.893058503327009E-010,&
        & 9.336586360508673E-013 /)
    integer,  dimension(nkpol_sgatie,nper_sgatie), parameter  :: mst_sgatie = (/ 0,0,0,0,0,0,0,1,1,1,1,1,1,2,2,2,2,2,3,3,3,3,4,4,4,&
                                                                                & 5,5,6,0,1,2,3,4,5,6,0,1,2,3,4,5,0,1,2,3,4,0,1,2,3,0,1,2,0,1,0/)
    !максимальная степень сжатия от расхода
    real, dimension(nkpol_sgatie_max), parameter :: kpol_sgatie_max = (/ -58.9598974870394,2.84819353954500,-5.793956209950129E-002,6.696787613820742E-004,&
        & -4.806629797945031E-006,2.192358965500257E-008,-6.201154973737457E-011,9.937347856322555E-014,-6.903738378205108E-017/)
    integer,  dimension(nkpol_sgatie_max), parameter  :: mst_sgatie_max = (/ 0,1,2,3,4,5,6,7,8 /)
    !минимальная степень сжатия от расхода
    real, dimension(nkpol_sgatie_min), parameter :: kpol_sgatie_min = (/271.398604947138,-17.2551149108950,0.477714817320597,&
        & -7.482113568823753E-003,7.248149037365788E-005,-4.444970301941715E-007,&
        & 1.684233500959385E-009,-3.602579439304682E-012,3.327621535570503E-015/)
    integer,  dimension(nkpol_sgatie_min), parameter  :: mst_sgatie_min =(/0,1,2,3,4,5,6,7,8/)
    !линия запирания: расход от частоты
    real, dimension(nkpol_zapiranie), parameter :: kpol_zapiranie = (/28635.1022677557, -198349.023433772, 569953.530072832, &
    & -866881.163253782, 737095.110663522, -332286.632732878, 62060.7897102170/)
    integer,  dimension(nkpol_zapiranie), parameter  :: mst_zapiranie = (/ 0, 1, 2, 3, 4, 5, 6/)
    !линия помпажа: расход от частоты
    real, dimension(nkpol_pompaj), parameter :: kpol_pompaj = (/ 6610.84954026143, -46225.3220794241, 134398.241888976, &
    & -206326.481525875, 176899.548159081, -80328.6101538380, 15096.3845711427/)
    integer,  dimension(nkpol_pompaj), parameter  :: mst_pompaj = (/ 0, 1, 2, 3, 4, 5, 6/)
    !к.п.д. от расхода и частоты
    real, dimension(nkpol_kpd), parameter :: kpol_kpd = (/ -183.612335385304, 1327.40712320045, -3904.78383497277, &
    & 6010.72569372929, -5109.40456219258, 1481.01434772030, 642.144492655562, -0.322329491043523, 1.16071835379725, &
    & -1.28632972185079, 0.139141053943114, 23.6387099897690, -32.2043232064208, 2.107547719005177E-003, -7.024369797135802E-003, &
    & 7.198432240718403E-003, -0.266166277363210, 0.387012218036991, -7.479889841080129E-006, 3.200377532657594E-005, 1.448063321633710E-003, &
    & -2.372002411580644E-003, -2.151628119975206E-008, -4.087589844544328E-006, 7.761846564416362E-006, 4.513210545352679E-009, &
    & -1.234109054673857E-008, 6.919903697017648E-012/)
    integer,  dimension(nkpol_kpd,nper_kpd), parameter  :: mst_kpd = (/ 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, &
    & 4, 4, 4, 5, 5, 6, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 3, 0, 1, 2, 0, 1, 0/)
    !линия запирания, степень сжатия от расхода
    real, dimension(nkpol_sgatie_min_zapiranie), parameter :: kpol_sgatie_min_zapiranie = (/-19.8870786153774,0.593103318477173,-6.968288942952134E-003,&
        & 4.354557178245967E-005,-1.519564966792630E-007,2.813521539685601E-010,-2.158016474199446E-013/)
    integer, dimension(nkpol_sgatie_min_zapiranie) :: mst_sgatie_min_zapiranie = (/0,1,2,3,4,5,6/)
    !линия помпажа, степень сжатия от расхода
    real,dimension(nkpol_sgatie_max_pompaj), parameter :: kpol_sgatie_max_pompaj = (/-102.276116470877,5.09535444187845,-0.102688857984385,&
        & 1.079526170524586E-003,-6.201965207225715E-006,1.831849675991496E-008,-2.140828530650738E-011/)
    integer, dimension(nkpol_sgatie_max_pompaj) :: mst_sgatie_max_pompaj = (/0,1,2,3,4,5,6/)
    !!!РАСЧЕТ!!!
    !от частоты находим расход по линии запирания
    zapiranie=polinom_one_per(kpol_zapiranie,mst_zapiranie,nkpol_zapiranie,freq)
    !от частоты находим расход по линии помпажа
    pompaj=polinom_one_per(kpol_pompaj,mst_pompaj,nkpol_pompaj,freq)
    q(1)=racxod-pompaj*kpompaj
    q(2)=zapiranie-racxod
    q(3)=pin-pmin
    if ((q(1) >= 0.).and.(q(2) >= 0.).and.(q(3) >= 0.)) then
        !от частоты и расхода находим степень сжатия
        call Polinom (sgatie,kpol_sgatie,mst_sgatie,nkpol_sgatie,nper_sgatie,racxod,freq)
        !от частоты и расхода находим к.п.д.
        call Polinom(kpd,kpol_kpd,mst_kpd,nkpol_kpd,nper_kpd,racxod,freq)
    else 
        sgatie=1.
        kpd=1.
    end if
    !от расхода находим максимальное сжатие
    !if ((87.5 < racxod).and.(racxod <= 137.518)) then
    !    sgatie_max=polinom_one_per(kpol_sgatie_max_pompaj,mst_sgatie_max_pompaj,nkpol_sgatie_max_pompaj,racxod)
    !else if ((137.518 < racxod).and.(racxod <= 250.256)) then
    !    sgatie_max=polinom_one_per(kpol_sgatie_max,mst_sgatie_max,nkpol_sgatie_max,racxod)
    !else
    !    sgatie_max=2.007
    !end if
    !!от расхода находим минимальное сжатие
    !if ((159. < racxod).and.(racxod <= 250.)) then
    !    sgatie_min=polinom_one_per(kpol_sgatie_min_zapiranie,mst_sgatie_min_zapiranie,nkpol_sgatie_min_zapiranie,racxod)
    !else if ((87.257 < racxod).and.(racxod <= 159.)) then
    !    sgatie_min=polinom_one_per(kpol_sgatie_min,mst_sgatie_min,nkpol_sgatie_min,racxod)
    !else
    !    sgatie_min=1.276
    !end if
    !Ограничения
    !write(*,'(A5,2x,A3,f16.8,2x,A6,f16.8,2x,A6,f16.8)')'16MBt','Sg=',sgatie,'SGMAX=',sgatie_max,'SGMIN=',sgatie_min
    !q(1)=sgatie_max-sgatie
    !q(2)=sgatie-sgatie_min

    q(4)=pmax-pin*sgatie

end subroutine gpa_2H637

subroutine gpa_NC63(q,kpd,sgatie,racxod,freq,kpompaj,pin)
!6,3 МВт НЦ-6,3В/76-2,2
    real, parameter :: pmin = 3.38329425, pmax = 7.453054
    integer,parameter :: nper_sgatie=2, nkpol_sgatie=28, nkpol_sgatie_max=6, nkpol_sgatie_min=7, nkpol_zapiranie=6, nkpol_pompaj=6
    integer,parameter :: nper_kpd=2, nkpol_kpd=28, nkpol_sgatie_min_zapiranie=7, nkpol_sgatie_max_pompaj=7
    real,intent(out) :: sgatie,kpd,q(4)
    real :: sgatie_max,sgatie_min,zapiranie,pompaj,polinom_one_per
    real,intent(in) :: racxod,freq,kpompaj,pin
    !степень сжатия от расхода и частоты
    real, dimension (nkpol_sgatie),parameter :: kpol_sgatie = (/ -68.3094419105049,560.289672070893,-1833.80365985733,3132.76612929410,&
        & -2954.20510032887,1496.72813699879,&
        & 132.218570767119,-1.13173182413583,6.50248016154857,-14.8860446745643,17.0018724336503,-11.7224317355505,&
        & -30.5759495326559,-6.732183709627869E-004,3.525405837917936E-003,&
        & -5.404466751524601E-003,5.357627193151574E-002,1.02765092216272,&
        & -3.135917143216720E-006,-1.654981559743951E-005,-5.864076217929862E-004,&
        & -1.703936618884658E-002,1.100669335962706E-007,3.639130026073703E-006,&
        & 1.570107025580745E-004,-9.430251216001653E-009,-7.642585590784897E-007,1.537132824767441E-009 /)
    integer,  dimension(nkpol_sgatie,nper_sgatie), parameter  :: mst_sgatie = (/0,0,0,0,0,0,0,1,1,1,1,1,1,2,2,2,2,2,&
        & 3,3,3,3,4,4,4,5,5,6,0,1,2,3,4,5,6,0,1,2,3,4,5,0,1,2,3,4,0,1,2,3,0,1,2,0,1,0 /)
    !максимальная степень сжатия от расхода
    real, dimension(nkpol_sgatie_max), parameter :: kpol_sgatie_max = (/5.87696149688408,-0.182825591150894,4.125124995369300E-003,&
        & -4.422878953073847E-005,2.213216551460056E-007,-4.297560102462568E-010 /)
    integer,  dimension(nkpol_sgatie_max), parameter  :: mst_sgatie_max = (/0,1,2,3,4,5/)
    !минимальная степень сжатия от расхода
    real, dimension(nkpol_sgatie_min), parameter :: kpol_sgatie_min = (/-3.90076273758454,0.645030392417193,-3.109016998637977E-002,&
        & 7.873773357852035E-004,-1.105744489787724E-005,8.148447335384848E-008,-2.465031595856968E-010/)
    integer,  dimension(nkpol_sgatie_min), parameter  :: mst_sgatie_min =(/0,1,2,3,4,5,6/)
    !линия запирания: расход от частоты
    real, dimension(nkpol_zapiranie), parameter :: kpol_zapiranie = (/-667.106411199048,4271.66647142000,-10203.9148358169,&
        & 12298.7813120814,-7299.69280336539,1708.92380620833/)
    integer,  dimension(nkpol_zapiranie), parameter  :: mst_zapiranie = (/0,1,2,3,4,5/)
    !линия помпажа: расход от частоты
    real, dimension(nkpol_pompaj), parameter :: kpol_pompaj = (/-272.304994661233,1690.41200883855,-3868.02474803565,&
        & 4540.60906365915,-2640.54562314642,608.613987770582/)
    integer,  dimension(nkpol_pompaj), parameter  :: mst_pompaj = (/0,1,2,3,4,5/)
    !к.п.д. от расхода и частоты
    real, dimension(nkpol_kpd), parameter :: kpol_kpd = (/1.66739325960506,-11.0355872311229,14.0460577780631,&
        & 26.6558066680898,-36.6700393853070,93.2789428135403,&
        & 98.8769095725969,4.646970861631997E-002,0.459158146173583,&
        & -2.23023313900704,1.53832739071524,-5.20189232373218,&
        & -8.47043262199721,-4.220500160280612E-003,1.877161878185723E-002,&
        & 8.605868417855203E-003,0.109014853876545,0.292706468061865,&
        & -1.737222974670025E-005,-2.789899506011237E-004,-1.376531679559919E-003,&
        & -5.200996037079045E-003,9.689170770183098E-007,9.778045754305189E-006,&
        & 5.121217441800617E-005,-2.734842674268582E-008,-2.678555077746858E-007,&
        & 5.809907932428399E-010/)
    integer,  dimension(nkpol_kpd,nper_kpd), parameter  :: mst_kpd = (/0,0,0,0,0,0,0,1,1,1,1,1,1,2,2,2,2,2,&
        & 3,3,3,3,4,4,4,5,5,6,0,1,2,3,4,5,6,0,1,2,3,4,5,0,1,2,3,4,0,1,2,3,0,1,2,0,1,0/)
    !линия запирания, степень сжатия от расхода
    real, dimension(nkpol_sgatie_min_zapiranie), parameter :: kpol_sgatie_min_zapiranie = (/-171.241935549217,10.3984852162512,-0.259256700311308,&
        & 3.417785664532135E-003,-2.509528931799605E-005,9.728389496896171E-008,-1.554073076975560E-010/)
    integer, dimension(nkpol_sgatie_min_zapiranie) :: mst_sgatie_min_zapiranie = (/0,1,2,3,4,5,6/)
    !линия помпажа, степень сжатия от расхода
    real,dimension(nkpol_sgatie_max_pompaj), parameter :: kpol_sgatie_max_pompaj = (/-331.125114903403,38.5521524441845,-1.85459292507367,&
        & 4.734955896484547E-002,-6.763305680312316E-004,5.126079421063940E-006,-1.610413032054716E-008/)
    integer, dimension(nkpol_sgatie_max_pompaj) :: mst_sgatie_max_pompaj = (/0,1,2,3,4,5,6/)
    !!!РАСЧЕТ!!!
    !от частоты находим расход по линии запирания
    zapiranie=polinom_one_per(kpol_zapiranie,mst_zapiranie,nkpol_zapiranie,freq)
    !от частоты находим расход по линии помпажа
    pompaj=polinom_one_per(kpol_pompaj,mst_pompaj,nkpol_pompaj,freq)
    q(1)=racxod-pompaj*kpompaj
    q(2)=zapiranie-racxod
    q(3)=pin-pmin
    if ((q(1) >= 0.).and.(q(2) >= 0.).and.(q(3) >= 0.)) then
        !от частоты и расхода находим степень сжатия
        call Polinom (sgatie,kpol_sgatie,mst_sgatie,nkpol_sgatie,nper_sgatie,racxod,freq)
        !от частоты и расхода находим к.п.д.
        call Polinom(kpd,kpol_kpd,mst_kpd,nkpol_kpd,nper_kpd,racxod,freq)
    else
        sgatie = 1.
        kpd=1.
    end if
    !!от расхода находим максимальное сжатие
    !if ((41.203 < racxod).and.(racxod <= 64.724)) then
    !    sgatie_max=polinom_one_per(kpol_sgatie_max_pompaj,mst_sgatie_max_pompaj,nkpol_sgatie_max_pompaj,racxod)
    !else if ((64.724 < racxod).and.(racxod <=119.308 )) then
    !    sgatie_max=polinom_one_per(kpol_sgatie_max,mst_sgatie_max,nkpol_sgatie_max,racxod)
    !else
    !    sgatie_max=2.124
    !end if
    !!от расхода находим минимальное сжатие
    !if (( 76.112< racxod).and.(racxod <=119.386 )) then
    !    sgatie_min=polinom_one_per(kpol_sgatie_min_zapiranie,mst_sgatie_min_zapiranie,nkpol_sgatie_min_zapiranie,racxod)
    !else if ((41.291< racxod).and.(racxod <=76.083 )) then
    !    sgatie_min=polinom_one_per(kpol_sgatie_min,mst_sgatie_min,nkpol_sgatie_min,racxod)
    !else
    !    sgatie_min=1.388
    !end if


    
    !Ограничения
    !q(1)=sgatie_max-sgatie
    !q(2)=sgatie-sgatie_min
    q(4)=pmax-pin*sgatie

end subroutine gpa_NC63    
    
end module raschet_gpa_new