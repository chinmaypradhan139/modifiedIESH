module marcus
implicit none
real*8, parameter :: pi=3.1416
real*8 :: gh,gama,KT,dG,Er,mass,omega
real*8 :: inpot(14)
real*8 rf,rb



contains
!.........................................
subroutine Marcus_Integration
!=============================================
! Integration of a function using Simpson rule 
!=============================================
implicit none
real*8 a, b, integral
integer n, i, kk,ip

open(25,file='fort.23')
do ip=1,14
read(25,*) inpot(ip)
enddo
close(25)

mass=2000
omega=2.0E-04
gh=inpot(8)
Er=mass*(omega**2)*(gh**2)/2
dG=inpot(9)
KT=inpot(10)
gama=inpot(4)

!write(*,*)gh,Er,dG,KT,gama

a = -0.01
b = 0.01

n = 100




do i=1,16
   call simpson(f,a,b,integral,n)
!   write (*,101) n, integral
  
   n = n*2
end do
rf=integral
!write(*,*) rf

n=100
do i=1,16
   call simpson(g,a,b,integral,n)
!   write (*,101) n, integral
   n = n*2
end do

rb=integral

!write(*,*) rb
write(17,*)((1/(rf+rb))/5000)*7 

!call Marcus_plot(rf,rb)

end subroutine

!........................................
  Function f(x)
!----------------------------------------
! Function for integration
!----------------------------------------
implicit none
double precision f, x


f=gama*(exp(x/KT))/(1+exp(x/KT))*exp(-(Er+dG+x)**2/(4*Er*KT))/sqrt(4*pi*Er*KT)



return
end

  Function g(x)
!----------------------------------------
! Function for integration
!----------------------------------------
implicit none
double precision g, x

g=gama*(1/(1+exp(x/KT)))*exp(-(Er-dG-x)**2/(4*Er*KT))/sqrt(4*pi*Er*KT)



return
end




 Subroutine simpson(f,a,b,integral,n)
!==========================================================
! Integration of f(x) on [a,b]
! Method: Simpson rule for n intervals  
! written by: Alex Godunov (October 2009)
!----------------------------------------------------------
! IN:
! f   - Function to integrate (supplied by a user)
! a	  - Lower limit of integration
! b	  - Upper limit of integration
! n   - number of intervals
! OUT:
! integral - Result of integration
!==========================================================
implicit none
double precision f, a, b, integral,s
double precision h, x
integer nint
integer n, i


! if n is odd we add +1 to make it even
if((n/2)*2.ne.n) n=n+1

! loop over n (number of intervals)
s = 0.0
h = (b-a)/dfloat(n)
do i=2, n-2, 2
   x   = a+dfloat(i)*h
   s = s + 2.0*f(x) + 4.0*f(x+h)
end do
integral = (s + f(a) + f(b) + 4.0*f(a+h))*h/3.0

return 

end subroutine simpson
!...................................................
subroutine Marcus_plot(kf,kb)
implicit none        
double precision time,population
double precision Kt
integer :: tim,omega_t, r
real*8, intent(in) :: kf,kb

!write(*,*) kf,kb
population=0
kt=kf+kb
!write(*,*) kb/kt,kf,kb
!write(17,*)((1/(kf+kb))/5000)*7 
omega_t=inpot(14)

do tim=0,5000*omega_t,100   
   time=real(tim)
   population=exp(-(kt*time-log(kf)))/kt+(kb/kt)

   write(45,*) time,population
  

   
enddo   


end subroutine
!..................................................

end module









