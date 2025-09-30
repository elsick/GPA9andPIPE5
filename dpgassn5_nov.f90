    subroutine dpgassn5(OGR,alphaneb,qnebal,Q63AL,Q72AL,Q82AL,Q102AL,OGR63AL,OGR72AL,OGR82AL,opg,kapvlsum,ntr,igg,itg,ipg,trpl,opg1,alpha,AL63,AL72,AL82,AL102,cher,pch,hom)
    implicit none
    interface
        subroutine dpgaspr1(ogr,gg,itg,ipg,opg1,trpl,dtr,cher)
            implicit none
            real :: ogr,gg,itg,ipg,opg1,trpl,dtr,cher
        end subroutine
    end interface
    
    real :: hom,pch,cher,trpl,ipg,itg,igg,kapvlsum,opg,alphaneb,dpsum,opgmin
    real, dimension(6) :: OGR,alpha,ntr,opg1,gg
    real, dimension(5) :: qnebal
    real, dimension(6), parameter :: dtr=(/0.63,0.72,0.82,1.02,1.22,1.42/)
    real, dimension(6), parameter :: kapvl1=(/7.2e07,8.1e07,9.9e07,1.53e08,2.16e08,2.43e08/)
	real, dimension(4) :: AL63,Q63AL
	real, dimension(3) :: AL72,Q72AL,OGR63AL
	real, dimension(2) :: AL82,Q82AL 
	real, dimension(2) :: OGR72AL 
    real :: Q102AL,AL102,OGR82AL 
    real :: qkoef
    integer :: i
    
    !    qkoef=1.
    
    qkoef=1.
    ntr=1.
    
	do i=1,4
		ntr(1)=ntr(1)+i*AL63(i)
	end do
	
    do i=1,4
		Q63AL(i)=AL63(i)*(1.-AL63(i))*qkoef
    end do
    
	do i=1,3
		OGR63AL(i)=AL63(i)-AL63(i+1)
    end do
    
	do i=1,3
		ntr(2)=ntr(2)+i*AL72(i)
	end do
	
	do i=1,3
		Q72AL(i)=AL72(i)*(1.-AL72(i))*qkoef
    end do
    
	do i=1,2
        OGR72AL(i)=AL72(i)-AL72(i+1)
    end do
	
    do i=1,2
		ntr(3)=ntr(3)+i*AL82(i)
	end do
	
	do i=1,2
		Q82AL(i)=AL82(i)*(1.-AL82(i))*qkoef
    end do
    
	OGR82AL=AL82(1)-AL82(2)
	
    ntr(4)=1.+AL102
	
    Q102AL=AL102*(1.-AL102)*qkoef

	do i=1,5
         QNEBAL(I)=ALPHA(I+1)*(1.-ALPHA(I+1))*qkoef
    end do

    kapvlsum=0.
    dpsum=0.
    gg=igg/ntr !‡ÒıÓ‰ Ì‡ 1 ÚÛ·Û
    cher=0.01 !¯ÂÓıÓ‚‡ÚÓÒÚ¸ ÚÛ·
    itg=300.
    
    do i=1,6
        call DPGASPR1(OGR(i),GG(i),itg,ipg,opg1(i),trpl,dtr(i),cher)
    end do
    
    
    kapvlsum=sum(trpl*kapvl1*ntr*alpha)
    dpsum=sum((ipg-opg1)*alpha)
    opg=ipg-dpsum
    alphaneb=1.-sum(alpha)
   
    
	    !write(*,*)'NTR:',ntr
    if (pch==3.) then
        30 Format(/5x, '=== ENTER PARAMETERS GASPIPE #',f3.0,' ===',//5x,'IGG:',f16.4,2x,'ITG:',f16.4,2x,'IPG:',f16.4,2x,'TRPL:',f16.4,&
           & //5x,'DTR:  ',6(f16.4),/5x,'ALPHA:',6(f16.4),/5x,'GG1:  ',6(f16.4),/5x,'OPG1: ',6(f16.4))
        40 format(/5x,'=== OUTPUT PARAMETERS GASPIPE #',f3.0,' ===',//5x,'ALPHANEB:',f16.4,2x,'OPG:',f16.4,2x,&
           & //5x,/5x,'OGR:',2x,6(f16.4),/5x,'QNEBAL:',5(f16.4),//5x,'KAPVLSUM:',f26.4,//5x,80('='))
        50 format(/5x,'AL63:',1x,4(f16.4),/5x,'AL72:',1x,3(f16.4),/5x,'AL82:',1x,2(f16.4),/5x,'AL102:',f16.4,//5x,&
           & 'Q63AL:',1x,4(f16.4),/5x,'Q72AL:',1x,3(f16.4),/5x,'Q82AL:',1x,2(f16.4),/5x,'Q102AL:',f16.4,//5x,&
           &'OGR63AL:',1x,3(f16.4),/5x,'OGR72AL:',1x,2(f16.4),/5x,'OGR82AL:',1x,1(f16.4),/5x,'Q102AL:',f16.4,//5x,'NTR:',6(f16.4))
        write(*,30)hom,igg,itg,ipg,trpl,dtr,alpha,gg,opg1
        write(*,40)hom,alphaneb,opg,OGR,qnebal,kapvlsum
        write(*,50)AL63,AL72,AL82,AL102,Q63AL,Q72AL,Q82AL,Q102AL,OGR63AL,OGR72AL,OGR82AL,q102al,ntr
    end if
    
    return
end subroutine dpgassn5
    
SUBROUTINE DPGASPR1(OGR,gg,itg,ipg,opg1,trpl,dtr,cher)
    implicit none
    
    interface
        subroutine reprg(re,pr,vm,zkm,zmum,gy,gm,ga,gu,gc,gw,gv,tg,pg,wg,d)
            implicit none
            real :: re,pr,vm,zkm,zmum,gy,gm,ga,gu,gc,gw,gv,tg,pg,wg,d
        end subroutine
        subroutine dptrp(DPTR,trpl,W,DTRP,PLO,RE,CHER,pech1)
            implicit none
            real :: cher
            real,intent(out) :: dptr
            real,intent(in) :: pech1,re,plo,trpl,w,dtrp
        end subroutine
    end interface 

    real :: OGR,gg,itg,ipg,opg1,trpl,dtr,cher
    REAL :: GY,GM,GA,GU,GC,GW
    real :: GV,PECH1,RE1,PR,VM,ZKM,ZMUM,WG,RE,PLO,DPYCH,dpy
    real :: p1,p2,DL,PCR,DLSUMM,DLTEC
    INTEGER I,jych

    jych=90  

      GY=0.0006*GG
      GM=0.938*GG
      GA=0.004*GG
      GU=0.000001*GG
      GC=0.000001*GG
      GW=0.000001*GG
      GV=0.000001*GG 
      
      pech1=0.

      dpy=(ipg-opg1)/jych
      p1=ipg
      p2=ipg-dpy
      DL=TRPL/jych
      DLSUMM=0.
      
      do i=1,jych
          
          PCR=P1-DPY/2
          call REPRG(RE1,PR,VM,ZKM,ZMUM,GY,GM,GA,GU,GC,GW,GV,itg,PCR,1.,dtr) 
          
          wg=vm*gm/(3.14*dtr**2/4.)
         
          call REPRG(RE,PR,VM,ZKM,ZMUM,GY,GM,GA,GU,GC,GW,GV,itg,PCR,wg,DTR)
          
          PLO=1/VM
          
          !if(pech1==2.) then
              !WRITE(*,*)'è´Æ‚≠Æ·‚Ï £†ß†=',PLO
              !WRITE(*,*)'ÇÔß™Æ·‚Ï £†ß†=',ZMUM
              !WRITE(*,*)'RE=',RE
          !end if
          
          CALL DPTRP(DPYCH,DL*1000.,WG,DTR,PLO,RE,CHER,pech1)
          
          DLTEC=DL*DPY/DPYCH
          
          !if (DLTEC>10000.) then
              !WRITE(*,*)'DLTEC,DL,DPY,DPYCH=',' ',DLTEC,DL,DPY,DPYCH
          !end if    
              
          DLSUMM=DLSUMM+DLTEC
          !write(*,*)'DDLSUM:',dlsumm,dltec
          P1=P2
          P2=P1-DPY
          !if(pech1==2.) WRITE(*,*)'DPYCH= ',DPYCH
      end do
      
      OGR=DLSUMM-TRPL
      
      return
END SUBROUTINE DPGASPR1
   
SUBROUTINE DPTRP(DPTR,trpl,W,DTRP,PLO,RE,CHER,pech1)
    real :: cher,cherot,QGIDRS
    real,intent(out) :: dptr
    real,intent(in) :: pech1,re,plo,trpl,w,dtrp
    
    CHER=0.01/1000.
    !cher=3e-005
    CHEROT=CHER/(DTRP)
    QGIDRS=0.067*(158./re+2*CHEROT)**0.2
    DPTR=(QGIDRS*(trpl/DTRP)*((PLO*W**2)/2.))/(0.98e005)
    if(pech1==2.) WRITE(*,*)'ÉàÑêÄÇãàóÖëäéÖ ëéèêéíàÇãÖçàÖ íêìÅéèêéÇéÑÄ- a‚a=',DPTR
    RETURN
    
END SUBROUTINE DPTRP
    
subroutine reprg(re,pr,vm,zkm,zmum,gy,gm,ga,gu,gc,gw,gv,tg,pg,wg,d)
    integer :: i,j
    real :: mum,km,kohc,cpm,cvm,rm1
    real, dimension(5) :: a1,a2,a3,a4,a5,a6,a7
    real, dimension(7) :: cp,cv,k,mu,m,y,v,tc,kz,kzz,sm,sf,sigma,epsk,y0
    real,dimension(5,7) :: a
    real, dimension(7,7) :: f
    real :: tt,tg,gcym,gy,gm,ga,gu,gc,gw,gv,sy,t,rm,pg,rom,vm,zkm,zmum,re,wg,d,pr
    data tc /304.2,190.7,126.2,133.0,373.6,647.27,33.3/
    data m /44.01,16.04,28.016,28.01,34.08,18.02,2.016/
    data sigma /3.91,3.758,3.798,3.69,3.623,2.641,2.827/
    data epsk /203.3,148.6,71.4,91.7,60.02,809.1,59.7/
    data a1 /3.8231419e-1,2.5207184e-4,-1.6633384e-7,7.6427112e-11,-2.0555466e-14/
    data a2 /3.7002972e-1,1.7125262e-4,5.8326777e-7,-1.1547931e-9,1.043852e-12/
    data a3 /3.0929091e-1,-5.3739164e-6,6.2620324e-8,-4.7710105e-11,1.543612e-14/
    data a4 /3.1004196e-1,9.3869493e-6,-3.1934747e-8,2.3601401e-10,-3.14685319e-13/
    data a5 /3.6002272e-1,5.1942602e-5,7.4235126e-8,-3.7660246e-11,-2.2581577e-14/
    data a6 /3.567226e-1,2.4795243e-5,5.7207222e-8,-3.5393369e-11,9.15388839e-15/
    data a7 /3.0502797e-1,3.8432404e-5,-9.7465055e-8,1.148019e-10,-4.37063309e-14/
     do i=1,5
        a(i,1)=a1(i)
    end do
    do i=1,5
        a(i,2)=a2(i)
    end do
    do i=1,5
        a(i,3)=a3(i)
    end do
    do i=1,5
        a(i,4)=a4(i)
    end do
    do i=1,5
        a(i,5)=a5(i)
    end do
    do i=1,5
        a(i,6)=a6(i)
    end do
    do i=1,5
        a(i,7)=a7(i)
    end do
    tt=tg
    gcym=gy+gm+ga+gu+gc+gw+gv
    y(1)=gy/gcym
    y(2)=gm/gcym
    y(3)=ga/gcym
    y(4)=gu/gcym
    y(5)=gc/gcym
    y(6)=gw/gcym
    y(7)=gv/gcym
    do i=1,7
        v(i)=22.41/m(i)
    end do
    sy=0
    do i=1,7
        sy=sy+y(i)/m(i)
    end do
    do i=1,7
        y0(i)=(y(i)/m(i))/sy
    end do
    !c  bõóàcãehàe teèãoemkoctà
    t=tt-273.15
    cpm=0
    cvm=0
    do i=1,7
        cp(i)=a(1,i)+a(2,i)*t+a(3,i)*t*t+a(4,i)*t*t*t+t*t*a(5,i)*t*t
        cp(i)=cp(i)*v(i)
    end do
    do i=1,7
       cv(i)=cp(i)-1.98725/m(i)
       cpm=y(i)*cp(i)+cpm
       cvm=cvm+y(i)*cv(i)
    end do
!c  bõóàcãehàe büákoctà komèohehtob
    kohc=0.000002669
    do i=1,7
        mu(i)=kohc*(m(i)*tt)**0.5/(sigma(i)**2/(0.697*(1+0.323*alog(tt/epsk(i)))))
    end do
!c  bõóàcãehàe teèãoèpoboÑhoctà
    do i=1,7
        k(i)=3600.*(1.3*cv(i)*m(i)+3.4-0.7*tc(i)/tt)*mu(i)/m(i)
    end do
    kz(2)=k(2)*(1./(1.+0.25*(cv(2)/(cp(2)-cv(2))-1.)))
    do 14 i=1,7
        if (i.eq.2) goto 151
        kz(i)=k(i)*(1./(1.+0.35*(cv(i)/(cp(i)-cv(i))-1.)))
        151 continue
    14 continue
    do i=1,7
        kzz(i)=k(i)-kz(i)
    end do
    do 16 i=1,7
        do 17 j=1,7
            if (i.eq.j) goto 161
            f(i,j)=(1+(mu(i)/mu(j))**0.5*(m(j)/m(i))**(0.25))**2/(8**0.5*(1+(m(i)/m(j)))**0.5)
            161 continue
        17 continue
    16 continue
    km=0.
    do 18 i=1,7
        sm(i)=0.
        sf(i)=0.
        do 19 j=1,7
            if (i.eq.j) goto 191
                sm(i)=sm(i)+(((m(i)+m(j))/(2*m(i)))**0.125)*f(i,j)*(y0(j)/y0(i))
                sf(i)=sf(i)+f(i,j)*(y0(j)/y0(i))
                191 continue
        19 continue
        km=km+kz(i)/(1+sm(i))+kzz(i)/(1+sf(i))
    18 continue
!c  bõóàcãehàe büákoctà cmecà Éaáob
     mum=0.
    do i=1,7
        mum=mum+mu(i)/(1+sf(i))
    end do
!c  bõóàcãehàe èãothoctà à yÑeãúhoÉo oÅöema
    rm=0.
    do i=1,7
        rm=rm+(y(i)*v(i))
    end do
    rm1=1./rm
!write(*,*)'rm1=tt=pg',rm1,tt,pg
    rom=273.0/tt*rm1*pg*(735./760.)
    vm=1./rom
    zkm=km
    zmum=mum
    re=wg*d*rom/mum
    pr=3600.*mum*cpm/km
    return
end subroutine reprg
