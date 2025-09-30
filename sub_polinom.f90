    subroutine Polinom (fn,kpol,mst,npol,nper,x1,x2)          
     integer,intent(in) :: npol,nper
     real,intent(in) :: x1,x2
     real, intent(out) :: fn
     integer, dimension (npol,nper),intent(in) :: mst
     real, dimension (nper) :: xn
     real, dimension (npol),intent(in) :: kpol
     real :: s
     integer :: i,j
     
     fn=0.
     xn(1)=x1
     xn(2)=x2
     do i=1,npol
        s=1
        do j=1,nper
            s=s*xn(j)**mst(i,j)
        end do
       fn=fn+s*kpol(i)
     end do
    end subroutine Polinom 
    
    real function polinom_one_per(kpol,mst,nkpol,x)
    !real :: polinom_one_per
    real, intent(in) :: x
    integer,intent(in) :: nkpol
    integer,dimension(nkpol),intent(in) :: mst
    real,dimension(nkpol),intent(in) :: kpol
    integer :: i
    polinom_one_per=0.
    do i=1,nkpol
        polinom_one_per=polinom_one_per+kpol(i)*x**mst(i)
    end do
    
    end function polinom_one_per