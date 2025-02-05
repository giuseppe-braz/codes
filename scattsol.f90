!This module contains the functions responsible for the criation of a selfdual solution (for the su2 theory)
module BPS
    implicit none

    contains
        

        !Main subroutine in this module, it creates (using the others functions in this module) a selfdual solution of the su2
        !theory, and gives back a vector containing the field in all of the discrete space.
        subroutine equation_su2(phi,x0,x1,dx,autodual)
            implicit none
            real(8), intent(inout), dimension(:,:) :: phi       !Field
            !real(8), intent(in), dimension(:,:) :: eta, inv_eta          
            !real(8), intent(in), dimension(:) :: params     
            real(8), intent(in) :: dx, x0, x1                   !Parameters, x0 gives the "position" of the kink, x1 is the border
                                                                !of the space, dx is the spacing of the lattice
            !real(8), dimension(:), allocatable :: k1, k2, k3, k4, dphi, aux2
            real(8) :: k1, k2, k3, k4, dphi, aux2
            real(8) :: aux
            integer, intent(in) :: autodual
            !integer, intent(in) :: tamanho, autodual            !autodual is a variable that tells the program if we want a kink or
                                                                !antikink (1 -> kink, -1 -> antikink)
            integer :: nx, nx_frente, nx_tras, sinal, i, j

            !Allocando os vetores do Runge-Kutta
            !allocate(k1(tamanho),k2(tamanho),k3(tamanho),k4(tamanho), dphi(tamanho))
            !allocate(aux2(tamanho))



            !Parametros do tamanho dos vetores        
            nx = (x1-x0)/dx                     !this gives the number of points in the lattice
            nx_frente = (x1-x0)/dx              !gives the number of points "forward" of the kink
            nx_tras = abs((-x1-x0)/dx)          !gives the number of points "back" of the kink   
                
            !if (tamanho.eq.1) then      

        !Now, we begin the runge-kutta method

            sinal = -1                              !here we will start constructing the solution backwards (from x0 -> -x1)

            do i = 1, (nx_tras)
              aux = x0 + (sinal*(i)*dx)             !current position
              !The next k's are numbers associated with the runge-kutta method
              k1 = dx*su2(phi(i,1),autodual)                 
              k2 = dx*su2(phi(i,1)+0.5d0*k1,autodual)
              k3 = dx*su2(phi(i,1)+0.5d0*k2,autodual)
              k4 = dx*su2(phi(i,1)+k3,autodual)
              dphi = (k1+2*k2+2*k3+k4)/6.d0
              phi(i+1,1) = phi(i,1)+sinal*dphi
            enddo
            
            
            call invert_array(phi,1,nx_tras)        !This puts the field vector in the right "direction", phi(i=0) <-> phi(x=-x1)


            sinal = 1                               !Now, constructing the solution forward 

            do i = nx_tras + 1, (nx_tras + nx_frente + 1)
              aux = x0 + (sinal*(i)*dx)
              k1 = dx*su2(phi(i,1),autodual)
              k2 = dx*su2(phi(i,1)+0.5d0*k1,autodual)
              k3 = dx*su2(phi(i,1)+0.5d0*k2,autodual)
              k4 = dx*su2(phi(i,1)+k3,autodual)
              dphi = (k1+2*k2+2*k3+k4)/6.d0
              phi(i+1,1) = phi(i,1)+sinal*dphi
            enddo


            !We have constructed the selfdual solution of the su2 theory centered in x0

            
        end subroutine equation_su2


   subroutine equation_su2_mod(phi,x0,x1,dx,autodual,b)
            implicit none
            real(8), intent(inout), dimension(:,:) :: phi       !Field
            !real(8), intent(in), dimension(:,:) :: eta, inv_eta          
            !real(8), intent(in), dimension(:) :: params     
            real(8), intent(in) :: dx, x0, x1, b                   !Parameters, x0 gives the "position" of the kink, x1 is the border
                                                                !of the space, dx is the spacing of the lattice
            !real(8), dimension(:), allocatable :: k1, k2, k3, k4, dphi, aux2
            real(8) :: k1, k2, k3, k4, dphi, aux2
            real(8) :: aux
            integer, intent(in) :: autodual
            !integer, intent(in) :: tamanho, autodual            !autodual is a variable that tells the program if we want a kink or
                                                                !antikink (1 -> kink, -1 -> antikink)
            integer :: nx, nx_frente, nx_tras, sinal, i, j

            !Allocando os vetores do Runge-Kutta
            !allocate(k1(tamanho),k2(tamanho),k3(tamanho),k4(tamanho), dphi(tamanho))
            !allocate(aux2(tamanho))


            !Parametros do tamanho dos vetores        
            nx = (x1-x0)/dx                     !this gives the number of points in the lattice
            nx_frente = (x1-x0)/dx              !gives the number of points "forward" of the kink
            nx_tras = abs((-x1-x0)/dx)          !gives the number of points "back" of the kink   
                
            !if (tamanho.eq.1) then      

        !Now, we begin the runge-kutta method

            sinal = -1                              !here we will start constructing the solution backwards (from x0 -> -x1)

            do i = 1, (nx_tras)
              aux = x0 + (sinal*(i)*dx)             !current position
              !The next k's are numbers associated with the runge-kutta method
              k1 = dx*su2_mod(phi(i,1),autodual,b)                 
              k2 = dx*su2_mod(phi(i,1)+0.5d0*k1,autodual,b)
              k3 = dx*su2_mod(phi(i,1)+0.5d0*k2,autodual,b)
              k4 = dx*su2_mod(phi(i,1)+k3,autodual,b)
              dphi = (k1+2*k2+2*k3+k4)/6.d0
              phi(i+1,1) = phi(i,1)+sinal*dphi
            enddo
            
            
            call invert_array(phi,1,nx_tras)        !This puts the field vector in the right "direction", phi(i=0) <-> phi(x=-x1)


            sinal = 1                               !Now, constructing the solution forward 

            do i = nx_tras + 1, (nx_tras + nx_frente + 1)
              aux = x0 + (sinal*(i)*dx)
              k1 = dx*su2_mod(phi(i,1),autodual,b)
              k2 = dx*su2_mod(phi(i,1)+0.5d0*k1,autodual,b)
              k3 = dx*su2_mod(phi(i,1)+0.5d0*k2,autodual,b)
              k4 = dx*su2_mod(phi(i,1)+k3,autodual,b)
              dphi = (k1+2*k2+2*k3+k4)/6.d0
              phi(i+1,1) = phi(i,1)+sinal*dphi
            enddo


            !We have constructed the selfdual solution of the su2 theory centered in x0

            
        end subroutine equation_su2_mod



         subroutine equation_su3(phi,x0,x1,dx,autodual)
            implicit none
            real(8), intent(inout), dimension(:,:) :: phi       !Field
            real(8), dimension(2,2) :: eta, inv_eta          
            real(8),  dimension(4) :: params     
            real(8), intent(in) :: dx, x0, x1                   !Parameters, x0 gives the "position" of the kink, x1 is the border
                                                                !of the space, dx is the spacing of the lattice
            real(8), dimension(2) :: k1, k2, k3, k4, dphi, aux2
            real(8) :: aux, det
            integer, intent(in) :: autodual
            !integer, intent(in) :: tamanho, autodual            !autodual is a variable that tells the program if we want a kink or
                                                                !antikink (1 -> kink, -1 -> antikink)
            integer :: nx, nx_frente, nx_tras, sinal, i, j

            !Allocando os vetores do Runge-Kutta
            !allocate(k1(tamanho),k2(tamanho),k3(tamanho),k4(tamanho), dphi(tamanho))

            params(1) = 0.1d0         !gamma1
            params(2) = 0.5d0         !gamma2  
            params(3) = 0.5d0       !gamma3
            params(4) = 0.5d0       !lambda

            !matrix
            eta(1,1) = 2.d0
            eta(2,2) = 2.d0
            eta(1,2) = -params(4)
            eta(2,1) = eta(1,2)

            !Inverse matrix
            det = eta(1,1)*eta(2,2) - eta(1,2)*eta(2,1)
            inv_eta(1,1) = eta(2,2)/det
            inv_eta(2,2) = eta(1,1)/det
            inv_eta(1,2) = (-1)*eta(1,2)/det
            inv_eta(2,1) = (-1)*eta(2,1)/det


            !Parametros do tamanho dos vetores        
            nx = (x1-x0)/dx                     !this gives the number of points in the lattice
            nx_frente = (x1-x0)/dx              !gives the number of points "forward" of the kink
            nx_tras = abs((-x1-x0)/dx)          !gives the number of points "back" of the kink   

            sinal = -1
            do i = 1, nx_tras
                aux = x0 + (sinal*i*dx)
                
                do j = 1, 2
                    aux2(j) = phi(i,j)
                enddo

                do j = 1, 2
                    k1(j) = su3(aux2,inv_eta,params,j)
                enddo
                do j = 1, 2
                    k2(j) = su3(aux2+(0.5d0*k1),inv_eta,params,j)
                enddo
                do j = 1, 2
                    k3(j) = su3(aux2+(0.5d0*k2),inv_eta,params,j)
                enddo
                do j = 1, 2
                    k4(j) = su3(aux2+(k3),inv_eta,params,j)
                enddo
                dphi = (k1+ 2*k2 +2*k3 +k4)/6.d0
                do j = 1, 2
                    phi(i+1,j) = phi(i,j) + sinal*dphi(j)
                enddo
            enddo

            call invert_array(phi,1,nx_tras)
            call invert_array(phi,2,nx_tras)

            sinal = 1

            do i = nx_tras, (nx_frente + nx_tras +1)
                aux = x0 + (sinal*i*dx)
                
                do j = 1, 2
                    aux2(j) = phi(i,j)
                enddo

                do j = 1, 2
                    k1(j) = su3(aux2,inv_eta,params,j)
                enddo
                do j = 1, 2
                    k2(j) = su3(aux2+(0.5d0*k1),inv_eta,params,j)
                enddo
                do j = 1, 2
                    k3(j) = su3(aux2+(0.5d0*k2),inv_eta,params,j)
                enddo
                do j = 1, 2
                    k4(j) = su3(aux2+(k3),inv_eta,params,j)
                enddo
                dphi = (k1+ 2*k2 +2*k3 +k4)/6.d0
                do j = 1, 2
                    phi(i+1,j) = phi(i,j) + sinal*dphi(j)
                enddo
            enddo


         end subroutine equation_su3




        !This function inverts the array
        subroutine invert_array(a,j,m)

          implicit none
          real(8), dimension(:,:), intent(inout) :: a
          integer :: head, tail, i, n, m, j
          real(8) :: aux

          n = m+1
          head = 1
          tail = n

          do while (head < tail)
            aux = a(head,j)
            a(head,j) = a(tail,j)
            a(tail,j) = aux
            head = head + 1
            tail = tail - 1
          end do

        end subroutine invert_array


        !This function is the "identity" of the theory, in this case is the su2
        real(8) function su2(phi,sinal) result(dphidt)
            implicit none
            real(8) :: phi
            integer :: sinal

            dphidt = 2.d0*dsin(0.5d0*phi)*sinal
        return
        end function

        real(8) function su2_mod(phi,sinal,b) result(dphidt)
            implicit none
            real(8) :: phi, b
            integer :: sinal

            dphidt = sinal*dsin(phi)*(1.d0 + (b*dcos(phi)))
        return
        end function
                  
        real(8) function su3(phi,inv_eta,params,n) result(dphidt)    
            real(8) :: lamb, gam1, gam2, gam3
            real(8), dimension(2,2) :: inv_eta
            real(8), dimension(4) :: params
            real(8), dimension(2) :: phi, U_phi
            integer :: i,j,n
                      !n define qual campo vamos avaliar a eq de 1 ordem
            

            gam1 = params(1)
            gam2 = params(2)
            gam3 = params(3)
            lamb = params(4)

            U_phi(1) = (-1)*gam1*dsin(phi(1)) - gam3*dsin(phi(1)-phi(2))    !Derivada do pre-potencial em relacao ao campo 1
            U_phi(2) = (-1)*gam2*dsin(phi(2)) + gam3*dsin(phi(1)-phi(2))    ! '' campo 2

            dphidt = inv_eta(n,1)*U_phi(1) + inv_eta(n,2)*U_phi(2)

        end function su3



end module



!This module is responsible of evolving the selfdual solution using the Verlet method, for the initial conditons we use the lorentz
!transformation method to relate the spatial derivative of the field in one reference frame with the temporal one in other.
module time_evol
    implicit none

    contains

        !main function, is the one resposible to evolve a initial solution
        subroutine evolution_su2(phi,dx,dt,Nx,Nt,c,gam,d)
            implicit none
            real(8), dimension(:,:), intent(inout) :: phi           !The field phi(i,j), i -> the time step (1,2,3), j -> spacial
                                                                    !lattice
            real(8), intent(in) :: dx, dt, c, gam                   !parameters, dx -> spacial lattice distance, dt -> temporal one
                                                                    !c -> velocity of the solution, gam -> lorentz gamma factor
            integer, intent(in) :: Nx, Nt, d                        !Nx -> number of points in the sapcial lattice, Nt -> number of
                                                                    !points in the temporal lattice, d -> parameter for the size of
                                                                    !the vector
            real(8), dimension(d) :: dxphi           !This array keeps the info about the derivatives of the field
            integer :: i, j, lim                                    
            real(8) :: x, dx_novo                                   

            !Agora vamos calcular a derivada (em x) do campo
            do i = 1, Nx
                dxphi(i) = D_x(phi,1,i,dx,Nx)             !We use the function D_x to determinate the first spacial derivative of
                                                            !the field
            enddo

                
            !We now "correct" the spacing of the spacial lattice because of the boost
            dx_novo = dx/gam

            lim = 0.5d0*(Nx-1)

            do i = 1, lim
                dxphi(i) = -c*gam*dxphi(i)      !Now dxphi becomes the temporal derivative of the moving solution (left kink)
            enddo

            do i = lim + 1, Nx
                dxphi(i) = c*gam*dxphi(i)       !Now dxphi becomes the remporal derivative of the moving solution (right kink)
            enddo


            !The first step of the time evolution is done using Euler-Cromer method
            do i = 1,Nx
                phi(2,i) = phi(1,i) + dt*dxphi(i)
            enddo

            open(2,file='teste.dat')        !Exporting to make a gif of the field
            open(3,file='teste4.dat')       !Exporting to make a gif of the energy of the solutions    

            !beginning the verlet method
            do i = 1, Nt
                do j = 1, Nx
                    x = (-0.5d0*(Nx-1) + (j-1))*dx_novo                                                     !Current position
                    phi(3,j) = 2*phi(2,j) - phi(1,j) + (dt*dt)*(D2_x(phi,2,j,dx_novo,Nx) - dsin(phi(2,j)))    !Updating the field in
                                                                                                            !all space
                    write(2,*) x, phi(2,j)                                                              !saving the field in the
                                                                                                        !current position
                    write(3,*) x, energ_su2(phi(2,j),D_x(phi,2,j,dx_novo,Nx),D_t(phi,2,j,dt))           !saving the energy in the
                                                                                                        !current position
                enddo
                write(2,*) ''
                write(2,*) ''
                write(3,*) ''
                write(3,*) ''
                do j = 1, Nx                !Saving the next field value for repeating the verlet method
                    phi(1,j) = phi(2,j)
                    phi(2,j) = phi(3,j)
                enddo
            enddo
            close(2)
            close(3)


        end subroutine


  subroutine evolution_su2_mod(phi,dx,dt,Nx,Nt,c,gam,d,b)
            implicit none
            real(8), dimension(:,:), intent(inout) :: phi           !The field phi(i,j), i -> the time step (1,2,3), j -> spacial
                                                                    !lattice
            real(8), intent(in) :: dx, dt, c, gam                   !parameters, dx -> spacial lattice distance, dt -> temporal one
                                                                    !c -> velocity of the solution, gam -> lorentz gamma factor
            integer, intent(in) :: Nx, Nt, d                        !Nx -> number of points in the sapcial lattice, Nt -> number of
                                                                    !points in the temporal lattice, d -> parameter for the size of
                                                                    !the vector
            real(8), dimension(d) :: dxphi           !This array keeps the info about the derivatives of the field
            integer :: i, j, lim                                    
            real(8) :: x, dx_novo                                   
            real(8) :: b

            !Agora vamos calcular a derivada (em x) do campo
            do i = 1, Nx
                dxphi(i) = D_x(phi,1,i,dx,Nx)             !We use the function D_x to determinate the first spacial derivative of
                                                            !the field
            enddo

                
            !We now "correct" the spacing of the spacial lattice because of the boost
            dx_novo = dx/gam

            lim = 0.5d0*(Nx-1)

            do i = 1, lim
                dxphi(i) = -c*gam*dxphi(i)      !Now dxphi becomes the temporal derivative of the moving solution (left kink)
            enddo

            do i = lim + 1, Nx
                dxphi(i) = c*gam*dxphi(i)       !Now dxphi becomes the remporal derivative of the moving solution (right kink)
            enddo


            !The first step of the time evolution is done using Euler-Cromer method
            do i = 1,Nx
                phi(2,i) = phi(1,i) + dt*dxphi(i)
            enddo

            open(2,file='teste.dat')        !Exporting to make a gif of the field
            open(3,file='teste4.dat')       !Exporting to make a gif of the energy of the solutions    

            !beginning the verlet method
            do i = 1, Nt
                do j = 1, Nx
                    x = (-0.5d0*(Nx-1) + (j-1))*dx_novo                                                     !Current position
                    phi(3,j) = 2*phi(2,j) - phi(1,j) + (dt*dt)*(D2_x(phi,2,j,dx_novo,Nx) - el_su2_mod(phi(2,j),b))    !Updating the field in
                                                                                                            !all space
                    write(2,*) x, phi(2,j)                                                              !saving the field in the
                                                                                                        !current position
                    write(3,*) x, energ_su2_mod(D_x(phi,2,j,dx_novo,Nx),D_t(phi,2,j,dt),pot_su2_mod(phi(2,j),b))           !saving the energy in the
                                                                                                        !current position
                enddo
                write(2,*) ''
                write(2,*) ''
                write(3,*) ''
                write(3,*) ''
                do j = 1, Nx                !Saving the next field value for repeating the verlet method
                    phi(1,j) = phi(2,j)
                    phi(2,j) = phi(3,j)
                enddo
            enddo
            close(2)
            close(3)


        end subroutine

        real(8) function el_su2_mod(phi,b)
            implicit none
            real(8) :: phi, b

            el_su2_mod = -b*dsin(phi) + 0.5d0*dsin(2*phi) + 3*b*dsin(3*phi) + 4*b*b*dsin(4*phi)
        return
        end

        real(8) function pot_su2_mod(phi,b)
            implicit none
            real(8) :: phi, b

            pot_su2_mod = 0.5d0*dsin(phi)*dsin(phi)*((1.d0 + 4*b*dcos(phi))**2)
        return
        end

        real(8) function energ_su2_mod(dxphi,dtphi,pot)
            implicit none
            real(8) :: dxphi, dtphi, pot

            energ_su2_mod = 0.5d0*(dxphi*dxphi + dtphi*dtphi) + pot

        return
        end




        subroutine evolution_su3(phi,dx,dt,Nx,Nt,c,gam,d,eta,inv_eta,tamanho,params)
            implicit none
            real(8), dimension(:,:), intent(inout) :: phi           !The field phi(i,j), i -> fields (1, ..., tamanho), j -> spacial
                                                                    !lattice
                                                                    
            real(8), dimension(:) :: params                                                        
            real(8), dimension(:,:) :: eta, inv_eta
            real(8), intent(in) :: dx, dt, c, gam                   !parameters, dx -> spacial lattice distance, dt -> temporal one
                                                                    !c -> velocity of the solution, gam -> lorentz gamma factor
            integer, intent(in) :: Nx, Nt, d, tamanho               !Nx -> number of points in the sapcial lattice, Nt -> number of
                                                                    !points in the temporal lattice, d -> parameter for the size of
                                                                    !the vector
            real(8), dimension(:,:), allocatable :: dxphi, phi1, phi2 
                                                                    
            integer :: i, j, lim, k, m                                    
            real(8) :: x, dx_novo                                   

            
            allocate(dxphi(tamanho,d), phi1(tamanho,d), phi2(tamanho,d))

            !Agora vamos calcular a derivada (em x) do campo
            do j = 1, tamanho
                do i = 1, Nx
                    dxphi(j,i) = D_x(phi,j,i,dx,Nx)             !We use the function D_x to determinate the first spacial derivative of
                enddo                                            !the field
            enddo

                
            !We now "correct" the spacing of the spacial lattice because of the boost
            dx_novo = dx/gam

            lim = 0.5d0*(Nx-1)

            do j = 1, tamanho
                do i = 1, lim
                    dxphi(j,i) = -c*gam*dxphi(j,i)      !Now dxphi becomes the temporal derivative of the moving solution (left kink)
                enddo
            enddo

            do j = 1, tamanho
                do i = lim + 1, Nx
                    dxphi(j,i) = c*gam*dxphi(j,i)       !Now dxphi becomes the remporal derivative of the moving solution (right kink)
                enddo
            enddo

            !The first step of the time evolution is done using Euler-Cromer method
            do j = 1, tamanho
                do i = 1,Nx
                    phi1(j,i) = phi(j,i) + dt*dxphi(j,i)
                enddo
            enddo

            open(2,file='teste.dat')        !Exporting to make a gif of the field
            open(3,file='teste4.dat')       !Exporting to make a gif of the energy of the solutions    

            !beginning the verlet method
            do i = 1, Nt
                do j = 1, Nx
                    x = (-0.5d0*(Nx-1) + (j-1))*dx_novo                                                     !Current position
                    do k = 1,tamanho
                        phi2(k,j) = 2*phi1(k,j) - phi(k,j) + (dt*dt)*(D2_x(phi1,k,j,dx_novo,Nx) - dsin(phi1(k,j)))    !Updating the field in
                    enddo                                                                                        !all space
                    write(2,*) x, (phi1(1,j), m=1,tamanho)                                                              !saving the field in the
                                                                                                        !current position
                    !write(3,*) x, energ_su2(phi(2,j),D_x(phi,2,j,dx_novo,Nx),D_t(phi,2,j,dt))           !saving the energy in the
                                                                                                        !current position
                enddo
                write(2,*) ''
                write(2,*) ''
                write(3,*) ''
                write(3,*) ''
                do j = 1, Nx                !Saving the next field value for repeating the verlet method
                    phi(1,j) = phi(2,j)
                    phi(2,j) = phi(3,j)
                enddo
            enddo
            close(2)
            close(3)


        end subroutine

        !Function that computes the first spacial derivative
        real(8) function D_x(phi,k,j,dx,Nx)
            implicit none
            real(8) :: dx
            real(8), intent(in), dimension(:,:) :: phi
            integer :: k, j, Nx

            if (j.eq.0) then
                D_x = (-137*phi(k,j)+300*phi(k,j+1)-300*phi(k,j+2)+200*phi(k,j+3)-75*phi(k,j+4)+12*phi(k,j+5))/(60*dx)
            else if (j.eq.Nx) then
                D_x = (137*phi(k,j)-300*phi(k,j-1)+300*phi(k,j-2)-200*phi(k,j-3)+75*phi(k,j-4)-12*phi(k,j-5))/(60*dx)
            else if((j.eq.1).or.(j.eq.Nx-1)) then
                D_x = (phi(k,j+1)-phi(k,j-1))/(2*dx)
            else
                D_x = (-phi(k,j+2)+8*phi(k,j+1)-8*phi(k,j-1)+phi(k,j-2))/(12*dx)
            endif

        return
        end

        !Function that computes the first temporal derivative
        real(8) function D_t(phi,k,j,dt)
            implicit none
            real(8) :: dt
            real(8), intent(in), dimension(:,:) :: phi
            integer :: k, j

            D_t = (phi(k+1,j)-phi(k-1,j))/(2*dt)

        return
        end

        !Second spatial derivative
        real(8) function D2_x(phi, k,  j, dx, Nx)
            implicit none
            real(8) :: dx
            real(8), intent(in), dimension(:,:) :: phi
            integer :: j, Nx, k

            if ((j.eq.2).or.(j.eq.(Nx-1))) then
                !D2_x = (15*3*phi(2,j) - 77*2*phi(2,j+1) + 107*2*phi(2,j+2) - 13*12*phi(2,j+3))/(12*(dx*dx)) + &
                !(61*phi(2,j+4)-10*phi(2,j+5))/(12.d0*(dx*dx)) - dsin(phi(2,j))
                D2_x = (phi(k,j-1) -2*phi(k,j) + phi(k,j+1))/(dx*dx)
                !D2_x = -dsin(phi(2,j))
            else if (j.eq.1) then
                D2_x = (45.d0*phi(k,j) - 154.d0*phi(k,j+1) + 214.d0*phi(k,j+2) - 156.d0*phi(k,j+3))/(12*dx*dx) + &
                    (61.d0*phi(k,j+4)-10.d0*phi(k,j+5))/(12.d0*dx*dx)
            else if (j.eq.Nx) then
                D2_x = (2.d0*phi(k,j) - 5.0d0*phi(k,j-1) + 4.d0*phi(k,j-2) - phi(k,j-3))/(dx*dx)
            else
                D2_x = (-phi(k,j-2) + 16.d0*phi(k,j-1) - 30.d0*phi(k,j) + 16.d0*phi(k,j+1) - phi(k,j+2))/(12*(dx*dx))
            endif

        return
        end

        !Energy of the su2 theory
        real(8) function energ_su2(phi,dxphi,dtphi)
            implicit none
            real(8) :: phi, dxphi, dtphi

            energ_su2 = 0.5d0*(dxphi*dxphi + dtphi*dtphi) + (1.d0-dcos(phi))

        return
        end

        real(8) function pot_su3(phi,inv_eta,params,j)
            implicit none
            real(8), dimension(:,:) :: phi
            real(8), dimension(:,:) :: inv_eta
            real(8), dimension(:) :: params
            integer, intent(in) :: j
            integer :: i, k
            real(8), dimension(2) :: dU
            real(8) :: v

            dU(1) = -params(1)*dsin(phi(1,j)) - params(3)*dsin(phi(1,j)-phi(2,j))
            dU(2) = -params(2)*dsin(phi(2,j)) + params(3)*dsin(phi(1,j)-phi(2,j))
            
            v = 0.d0

            do i = 1, 2
                do k = 1,2
                    v = v + 0.5d0*inv_eta(i,k)*dU(k)
                enddo
            enddo

            pot_su3 = v

        return
        end



end module

!there is important comments in the modules of the program, be sure to read them

!now the main program, this only has the function of giving some changing parameters and calling the important subroutines
program main

    use BPS
    use time_evol
    implicit none

    real(8), parameter :: pi = dacos(-1.d0)
    integer, parameter :: d = 200000            !size of the vectors

    real(8), dimension(3,d) :: phi, dxphi       
    real(8), dimension(:,:), allocatable :: phi0, eta, inv_eta      !phi0 is the initial position of the kink
    real(8), dimension(4) :: params             
    real(8) :: dx, dt, x, t                     
    real(8) :: alpha, c, gam, aux, L, Tmax, x0, det
    integer :: i, j, Nx, Nt, tamanho, autodual
    real(8), dimension(7) :: b


    !Decidindo o caso (su2, su3)
    tamanho = 1
    allocate(phi0(d,tamanho), eta(tamanho,tamanho), inv_eta(tamanho,tamanho))
    

    !Espaçamento da rede (t,x)
    dx = 0.1d0
    dt = dx/2.d0
    alpha = dt/dx

    !Tamanho da Rede (t,x)
    L = 50
    Tmax = 100

    !Número de entradas dos vetores
    Nx = (2*L)/dx +1        !Numero de pontos espaciais da rede
    Nt = ceiling(Tmax/dt)   !Numero de pontos temporais da rede

    !Parâmetros para a condição inicial do campo
    c = 0.2d0
    gam = 1/dsqrt(1 - c*c)


    !Dados iniciais via eq de bps para criar a solucao da esquerda
    autodual = 1            !Decide se kink ou antikink

    if (tamanho.eq.1) then
        !phi0(1,1) = pi + 0.01d0
        phi0(1,1) = 0.5d0*pi + 0.001d0
    else if (tamanho.eq.2) then
        phi0(1,1) = 0.1d0
        phi0(1,2) = 3.1095d0
    endif

    x0 = -35.d0             !posicao inicial do bixo

    !call equation(phi0,tamanho,eta,inv_eta,params,x0,L,dx,autodual)

    !b(1) = -2.d0
    !b(2) = -1.d0
    !b(3) = -0.5d0
    !b(4) = 0.d0
    !b(5) = 0.5d0
    !b(6) = 1.d0
    !b(7) = 2.d0
    
    !open(1,file='trip1.dat')
    !open(2,file='trip2.dat')
    !open(3,file='trip3.dat')
    !open(4,file='trip4.dat')
    !open(5,file='trip5.dat')
    !open(6,file='trip6.dat')
    !open(7,file='trip7.dat')


    !do i = 1,7
    !    phi0 = 0.d0
    !    phi0(1,1) = 0.5d0*pi + 0.01d0
    !    call equation_su2_mod(phi0,x0,L,dx,autodual,b(i))
    !    do j = 1, Nx
    !        write(i,*) (-L + j*dx), phi0(j,1)
    !    enddo
    !enddo

    !close(1)
    !close(2)
    !close(3)
    !close(4)
    !close(5)
    !close(6)
    !close(7)
        
    call equation_su2_mod(phi0,x0,L,dx,autodual,-2.d0)
    !call equation_su2(phi0,x0,L,dx,autodual)

    do i = 1, Nx
        phi(1,i) = phi0(i,1)
    enddo

    !open(2,file='teste.dat')
    !do i = 1, Nx
    !    write(2,*) (-L + i*dx), phi0(i,1), 4.d0*datan(dexp(-L + i*dx))
    !enddo

    !Dados inicias para criar a solucao da direita
    autodual = 1
    phi0 = 0.d0
    if (tamanho.eq.1) then
        !phi0(1,1) = pi + 0.1d0
        phi0(1,1) = 0.5d0*pi + 0.001d0
    else if (tamanho.eq.2) then
        phi0(1,1) = pi + 0.1d0
        phi0(1,2) = 1.3d0
    endif

    x0 = 35.d0

    call equation_su2_mod(phi0,x0,L,dx,autodual,-2.d0)

    do i = 1, Nx
        phi(1,i) = phi(1,i) + phi0(i,1)
    enddo

    !call evolution_su2_mod(phi,dx,dt,Nx,Nt,c,gam,d,-2.d0) 

    open(2,file='teste2.dat')
    do i = 1,Nx
        write(2,*) (-L + i*dx), phi(1,i)
    enddo
    close(2)

    !call evolution_su2(phi,dx,dt,Nx,Nt,c,gam,d)

end program main
