subroutine gas_propetries(qtek,qrnsumm,rosm,rsm,z,sost_dol,p,t,igg)
    !igg - массовый расход газа на входе, кг/с
    !g - расходы компонент газа в кмоль/с
    !GB - расход компонент газа на входе, кг/с
    !m - молярная масса компоненты
    !y - мольные доли компонент
    !tkr,pkr - критические температура и давление компонентов газа
    !rsm - газовая постоянная газа
    !rosm - плотность газа кг/м3
    !z - коэффициент сжимаемости газа
    !qtek - объемный расход, м3/мин
    !muu- динамическая вязкость газа, Па*с
    real, dimension(7), parameter :: m=(/16.043,30.07,44.097,58.123,72.15,44.01,28.0135/)!молекулярные массы компонентов газа CH4,C2H6,C3H8,C4H10,C5H12,CO2,N2
    real, dimension(7), parameter :: rost=(/0.6682,1.2601,1.8641,2.4956,3.228,1.8393,1.1649/)!плотность компонентов газа при стандартных условиях CH4,C2H6,C3H8,C4H10,C5H12,CO2,N2, (кг/м3)
    real, dimension(7), parameter :: tkr=(/190.66,305.46,369.9,425.2,469.5,304.26,126.2/)!критическая температура компонентов газа,CH4,C2H6,C3H8,C4H10,C5H12,CO2,N2, (К)
    real, dimension(7), parameter :: pkr=(/4.64,4.884,4.255,3.799,3.373,7.386,3.394/)!критическое давление компонентов газа, CH4,C2H6,C3H8,C4H10,C5H12,CO2,N2, (МПа)
    real, dimension(7), parameter :: qrn=(/35850.,63850.,91300.,118700.,146200.,0.,0./)
    real, dimension(7) ::g,y,GB
    real, intent(in) :: igg,p,t
    real, dimension(7), intent(in) :: sost_dol
    real :: gsumm,pkrsumm,tkrsumm,tau,msr,ppr,tpr,muu0,muu,B1,B2,B3
    real, intent(out) :: qtek,rsm,rosm,z,qrnsumm
    real,parameter :: r0=8314.3
    
    GB=sost_dol*igg*1000.
    g=GB/m
    gsumm=sum(g)
    y=g/gsumm
    pkrsumm=sum(pkr*y)
    tkrsumm=sum(tkr*y)
    ppr=p/pkrsumm
    tpr=t/tkrsumm
    tau=1.-1.68*tpr+0.78*tpr**2.+0.0107*tpr**3.
    z=1-0.0241*ppr/tau
    msr=sum(m*y)
    rsm=r0/msr
    rosm=p*1000000./(z*rsm*t)
    qtek=gsumm/rosm
    qrnsumm=sum(qrn*sost_dol)
    B1=-0.67+2.36/tpr-1.93/tpr**2.
    B2=0.8-2.89/tpr+2.65/tpr**2.
    B3=-0.1+0.354/tpr-0.314/tpr**2.
    muu0=(1.81+5.95*tpr)*1.e-06
    muu=muu0*(1+B1*ppr+B2*ppr**2.+B3*ppr**3.)
end subroutine gas_propetries
    
