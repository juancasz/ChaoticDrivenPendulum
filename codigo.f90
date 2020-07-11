!Programa para resolver ecuación por RK4
! x=theta: Posición angular
! x'=theta'=w: Velocidad angular
! x''=theta'' : Aceleración angular
! x''(t)+qx'(t)+sin(x(t))=bcos(wo*t)     Ecuación diferencial no lineal péndulo simple caos forzado
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
    real :: q,b !constantes ecuación diferencial
    real :: h !paso
    real :: T !Periodo
    real,dimension(3) :: y_inf !valor inicial Y(ti)
    real,dimension(3) :: y_sup !siguiente valor por RK4 Y(ti+h)
    integer :: i,j !Contadores
    real :: t_ini,t_end !Rango  de tiempo en el que se desea evaluar solución
    integer :: n !Número de puntos
    real,parameter :: PI=4.*ATAN(1.) !PI=3.14159...
    real,allocatable :: array(:,:)  !Arreglo para guardar datos y enviar a bloc de notas
    real,allocatable :: poincare(:,:) !Arreglo para guardar datos sección Poincaré
    integer :: status !bandera

    j=0 !Iniciar contador datos Poincaré

    wo=0.6667 !Frecuencia angular
    q=2
    b=1.15

    T=(2*PI)/wo  !Periodo
    h=T/50.  !Paso 1/50 periodo

    !Resolver para intervalo entre 0 y 100T
    t_ini=0.
    t_end=100.*T

    !VALOR INICIAL  Y(0) = (0 )
    !                      -PI/2
    !                       0

    y_inf = (/0.,-PI/2,0./)

    n=INT((t_end-t_ini)/h) !Número de datos

    allocate(array(n,3),STAT=status)
    allocate(poincare(100,3),STAT=status)

!***************************************************
    !Tabla imprimiendo t, x y w
    !write(*,*) "ECUACIÓN DIFERENCIAL NO LINEAL PÉNDULO SIMPLE FORZADO"
    ! write(*,1)"t","x(t)","x'(t)"
    !1   format(3X,A1,7X,A4,7X,A5)

   ! write(*,2) y_inf  !Valor inicial
!********************************************************

    array(1,:)=y_inf  !Condición inicial: primeros datos solución

    do i=2,n
        call RK4(y_inf,wo,q,b,y_sup)

        y_inf = y_sup

        !theta > PI o  theta<PI  se retorna a [-PI,PI]
        if(y_sup(2) > PI) then
         y_sup(2)=y_sup(2)-AINT((y_sup(2)+PI)/(2.*PI))*(2*PI)
        else if (y_sup(2)< -PI) then
         y_sup(2)=-ABS(y_sup(2))+AINT((ABS(y_sup(2))+PI)/(2*PI))*(2*PI)
        end if

        !Cada que se completa un periodo se guarda el dato para graficar sección de Poincaré
        if (MOD(i-1,50)==0.) then
           j=j+1
           poincare(j,:)=y_sup
        end if

        array(i,:)=y_sup


      ! write(*,2) y_sup   En caso de que se quiera imprimir en consola
    end do



    OPEN(unit=1,file="data.txt") !Guardar todos los datos en data. txt
        DO i=1,n
          WRITE(1,2) (array(i,j),j=1,3)
        END DO
    CLOSE(1)

    OPEN(unit=2,file="poincare.txt") !Guardar datos gráfico de Poincaré en poincare.txt
        DO i=1,99
          WRITE(2,2) (poincare(i,j),j=1,3)
        END DO
    CLOSE(1)

2 format(X,F10.2,X,F10.3,X,F10.3)
    write(*,*) "Datos leídos", n
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
   real, dimension(3) :: F1,F2,F3,F4 !Términos RungeKutta

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

!Transformación F
subroutine F(y_in,wo,q,b,y_out)
  implicit none
  real,intent(in) :: wo,q,b
  real,intent(in),dimension(3) :: y_in
  real,intent(out),dimension(3)::y_out

  y_out(1)=1
  y_out(2)=y_in(3)
  y_out(3)=-q*y_in(3)-SIN(y_in(2))+b*COS(wo*y_in(1))

end subroutine
