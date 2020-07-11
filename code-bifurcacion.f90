!Programa para resolver ecuaci�n por RK4
! x=theta: Posici�n angular
! x'=theta'=w: Velocidad angular
! x''=theta'' : Aceleraci�n angular
! x''(t)+qx'(t)+sin(x(t))=bcos(wo*t)     Ecuaci�n diferencial no lineal p�ndulo simple caos forzado
! q=k/m
! g/l = 1
! b= f0/ml
! x(0)=-PI/2
! x'(0)=w(0)
! Cambio de variable y1=t, y2=x, y3=x'
! Se tiene sistema
!             y1'               y1              1
!   Y'(t) = ( y2') = F(Y(t))= F(y2  ) = (       y3              )
!             y3'               y3      -qy3-sin(y2)+bcos(wo*y1)
program pendulo_caos
    implicit none
    real :: wo !frecuencia angular
    real :: q,b !constantes ecuaci�n diferencial
    real :: h !paso
    real :: T !Periodo
    real,dimension(3) :: y_inf !valor inicial Y(ti)
    real,dimension(3) :: y_sup !siguiente valor por RK4 Y(ti+h)
    real,dimension(4)::by_sup !A los tres datos de t,theta y w se les adjunta la constante b
    integer :: i,j,k !Contadores
    real :: t_ini,t_end !Rango  de tiempo en el que se desea evaluar soluci�n
    integer :: s !N�mero de periodos dela fuerza externa
    real,parameter :: PI=4.*ATAN(1.) !PI=3.14159...
    real,allocatable :: poincare(:,:) !Arreglo para guardar datos secci�n Poincar�
    real :: b_end, b_ini,b_step !Realizar barrido par�metro b
    integer :: m !N�mero de datos b a usar
    integer :: x !Tama�o matriz datos Poincar�
    integer :: n !N�merto total de muestreos
    integer :: l
    integer :: status !bandera

    k=0 !Iniciar contador datos Poincar�
    l=0

    s=60 !60 Periodos

    wo=2./3. !Frecuencia angular
    q=0.5



    T=(2*PI)/wo  !Periodo
    h=T/50.  !Paso 1/50 periodo

    b_ini=1.0
    b_end=1.501
    b_step=0.001

    m=INT((b_end-b_ini)/b_step)

    x=(m+1)*(s-10) !Se eliminan los primeros 10 puntos de Poincare

    !Resolver para intervalo entre 0 y 100T
    t_ini=0.
    t_end=s*T

    !VALOR INICIAL  Y(0) = (0 )
    !                      -PI/2
    !                       0

    y_inf = (/0.,-PI/2,0./)

    n=INT((t_end-t_ini)/h) !Muestreo

    allocate(poincare(x,4),STAT=status)


    do i=1,m+1
        b = b_ini+(i-1)*b_step
        do j=1,n
            call RK4(y_inf,wo,q,b,y_sup)
            y_inf = y_sup
            !Cada que se completa un periodo se guarda el dato para graficar secci�n de Poincar�
            if (MOD(j,50)==0.) then

               k=k+1

               if(k>10) then     !Se eliminan los primeros 10 puntos de Poincar�

                   !theta > PI o  theta<PI  se retorna a [-PI,PI]
                    if(y_sup(2) > PI) then
                        y_sup(2)=y_sup(2)-AINT((y_sup(2)+PI)/(2.*PI))*(2*PI)
                    else if (y_sup(2)< -PI) then
                        y_sup(2)=-ABS(y_sup(2))+AINT((ABS(y_sup(2))+PI)/(2*PI))*(2*PI)
                    end if

                   l=l+1

                   by_sup(1)=b
                   by_sup(2)=y_sup(1) !t
                   by_sup(3)=y_sup(2) !theta
                   by_sup(4)=y_sup(3) !w
                   poincare(l,:)=by_sup
               end if

           end if

        end do

        k=0
        y_inf = (/0.,-PI/2,0./)

    end do



    OPEN(unit=1,file="poincare.txt") !Guardar datos gr�fico de Poincar� en poincare.txt
        DO i=1,x
          WRITE(1,2) (poincare(i,j),j=1,4)
        END DO
    CLOSE(1)

2 format(X,F10.3,X,F10.3,X,F10.3,X,F10.3)
   write(*,*) "Datos le�dos", x
end program

!RungeKutta de orden 4
subroutine RK4(y_inf,wo,q,b,y_sup)
   implicit none
   real,intent(in) :: wo,q,b
   real, intent(in), dimension(3) :: y_inf !Y(ti)
   real, intent(out),dimension(3) :: y_sup !Y(ti+h)
   real :: T !Periodo
   real :: h !Paso
   real,parameter :: PI=4.*ATAN(1.) !PI=3.14159...
   real, dimension(3) :: y_out ! y_out=F(y_in)
   real, dimension(3) :: F1,F2,F3,F4 !T�rminos RungeKutta

   T=(2*PI)/wo  !Periodo
   h=T/50.  !Paso 1/50 periodo

   call F(y_inf,wo,q,b,y_out)
   F1 = h*y_out
   call F(y_inf+0.5*F1,wo,q,b,y_out)
   F2 = h*y_out
   call F(y_inf+0.5*F2,wo,q,b,y_out)
   F3 = h*y_out
   call F(y_inf+F3,wo,q,b,y_out)
   F4 = h*y_out

   y_sup = y_inf + (1./6.)*(F1+2*F2+2*F3+F4)

end subroutine

!Transformaci�n F
subroutine F(y_in,wo,q,b,y_out)
  implicit none
  real,intent(in) :: wo,q,b
  real,intent(in),dimension(3) :: y_in
  real,intent(out),dimension(3)::y_out

  y_out(1)=1
  y_out(2)=y_in(3)
  y_out(3)=-q*y_in(3)-SIN(y_in(2))+b*COS(wo*y_in(1))

end subroutine
