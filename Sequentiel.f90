Program Projet
  !use fonction
  !use Gradient
  Implicit NONE


  ! PARAMETRE DU SYSTEME
  Integer:: Nx, Ny, i, Maxiter, j, kt, k, n, l, Nl
  Real*8:: Lx, Ly, dx, dy, dt, tfinal, D, t
  Real*8,dimension(:),allocatable:: C_0, C, b, kappa, r, d1, W
  Real*8,dimension(:,:),allocatable::A
  Real*8:: coeff_a, coeff_b, coeff_c
  real:: t1, t2

  REAL(KIND=8):: epsilon, alpha, beta, residu, drl, dwl, betal, residul

  call CPU_TIME( t1 )


  ! INITIALISATION

  OPEN(10,file='data', form='formatted', status='old')
  READ(10,*)
  READ(10,*) Lx, Ly
  READ(10,*)
  READ(10,*) D
  READ(10,*)
  READ(10,*) Nx, Ny
  READ(10,*)
  READ(10,*) tfinal
  CLOSE(10)

  Allocate(C_0(1:Nx*Ny), C(1:Nx*Ny), b(1:Nx*Ny), W(1:Nx*Ny), kappa(1:Nx*Ny), r(1:Nx*Ny), d1(1:Nx*Ny))
  Allocate(A(1:Nx*Ny,1:Ny*Nx)) !expliquer d'ou vient la taille des vecteurs

  dx=Lx/(Nx+1)
  dy=Ly/(Ny+1)
  dt=0.9d0*dx*dx/(4.d0*D)
  !dt=1.d0
  coeff_c=(1.0d0 - 2.0d0*D*dt/(dy*dy) - 2.0d0*D*dt/(dx*dx))
  !coeff_a=2.0d0*D*dt/(dy*dy) + 2.0d0*D*dt/(dx*dx)
  coeff_a= (D*dt/(dx*dx))
  coeff_b= (D*dt/(dy*dy))

  print*, coeff_a, coeff_b, coeff_c

  n = Nx*Ny

  !!!!!!!!!!!!!!!!!!!!!!!!!
  !   Contruct de A       !
  !!!!!!!!!!!!!!!!!!!!!!!!!

  A=0.d0

  A(1,1)=coeff_c
  A(1,2)=coeff_a
  A(2,1)=coeff_a
  A(n,n)=coeff_c
  A(n,n-1)=coeff_a
  A(n-1,n)=coeff_a

  do i=2,n-1
    do j=2,n-1
      if (i==j) then
        A(i,i)=coeff_c
        A(i,i-1)=coeff_a
        A(i,i+1)=coeff_a
      end if
    end do
  end do

  do j=0,n-Nx
    A(1+j,Ny+1+j)=coeff_b
  end do

  do i=0,n-Ny
    A(Nx+1+i,1+i)=coeff_b
  end do

  ! do i=1, n
  !   print *, A(i,:)
  ! end do

  Maxiter = int(tfinal/dt)

  !On donne la condition initiale à C
  C(:)= 1.0d0 ! <--- A FAIRE


  !On définit b
  b=0.d0 ! <--- A FAIRE


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                    Gradient Conjugué
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !initialisation Gradient conjugue

  W=matmul(A,C)

  DO i=1,N
    kappa(i) = C(i)
    r(i)     = W(i) - b(i)
    d1(i)    = W(i) - b(i)
  END DO

  residu = 0.d0
  DO i=1,N
    residu = residu + r(i)*r(i)
  ENDDO

  ! boucle du Gradient conjugue

  l=1
  DO WHILE ((l<10000).AND.( SQRT(residu) .ge. 0.0001 ))

    DO i=1,N
      W=matmul(A,d1)
    ENDDO
    drl = 0.0d0
    dwl = 0.0d0
    DO i = 1, N
      drl = drl + d1(i)*r(i)
      dwl = dwl + d1(i)*w(i)
    ENDDO
    alpha = drl/dwl
    DO i=1,N
      kappa(i) = kappa(i) - alpha*d1(i)
      r(i) = r(i) - alpha*W(i)
    END DO
    betal=0.d0

    DO i=1,N
      betal=betal+ r(i)*r(i)
    ENDDO
    beta = betal/residu
    DO i=1,N
      d1(i) = r(i) + beta*d1(i)
    ENDDO
    residu = 0.d0
    DO i=1,N
      residu = residu+r(i)*r(i)
    ENDDO
    l=l+1

    !Fin Gradient conjugue

  ENDDO
  print*,'l',l,'residu',residu
  DO i=1,N
    C(i)=kappa(i)
  ENDDO
  t = t + dt;

  ! Fin boucle en temps

  k=1
  DO j=0,Ny
    !do i=1,Nx
    print*, C(j:j+Nx)
    k=k+1
    !END DO
  end do

  Deallocate(C_0,C,b)
  call CPU_TIME( t2 )
  print *,t2 - t1

end program
