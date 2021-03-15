module automatocelular
    !implicit none
    
    contains

        !função de geração, com distribuição de n focos ao longo da simulação
        function geracao_inicial(m) result(matrizgerada)
            use OMP_LIB
            implicit none
            
            !tamanho da população
            INteger, INTENT(IN) :: m
            integer :: i
            integer, parameter :: n = 25
            !matriz final e posições dos focos
            integer :: matrizgerada(m,m), posicoes(2,n)
            real :: r(2,n)

            !posições são aleatórias e transoformadas em coordenadas
            call RANDOM_NUMBER(r)
            posicoes = r*m
            
            !colocação dos focos
            !$OMP PARALLEL DO SIMD
            do i = 1,n
                matrizgerada(posicoes(1,i),posicoes(2,i)) = 1
            end do
            !$OMP END PARALLEL DO SIMD

        end function geracao_inicial

        !geração do timer que ativa a infecção
        function geracao_timer(matrizinit,m) result(matriztimer)
            USE OMP_LIB
            implicit none

            integer, intent(in) :: m
            integer, intent(in) :: matrizinit(m,m)
            integer :: matriztimer(m,m)
            real :: r(m,m)

            !geração do timer, com mais valores próximos de 0
            call RANDOM_NUMBER(r)
            matriztimer = matrizinit*(r**2)*60

        end function geracao_timer

        !atualização do timer
        function atualizacao_timer(matrizinit,matriztimerinit,m) result(matriztimer)
            USE OMP_LIB
            implicit none

            integer, intent(in) :: m
            integer :: i,j
            integer, intent(in) :: matrizinit(m,m)
            integer, intent(in) :: matriztimerinit(m,m)
            integer :: matriztimer(m,m)

            matriztimer = matriztimerinit 

            !atualização propriamente dita
            !$OMP PARALLEL DO SIMD
            do i = 1,m
                do j = 1,m
                    if ( matrizinit(j,i) >= 1 ) then
                        matriztimer(j,i) = matriztimer(j,i) + 1
                    end if
                end do
            end do
            !$OMP END PARALLEL DO SIMD

        end function atualizacao_timer

        !parte legal, infecção
        function calculo_infeccao(matrix,matriztimer,m) result(matrizfim)
            USE OMP_LIB
            implicit none

            !tamanho, matriz de população e timer
            integer, intent(in) :: m
            integer, intent(in) :: matrix(m,m)
            integer, intent(in) :: matriztimer(m,m)
            !matriz para transformação em um espaço toroidal
            integer :: matrizt_toroide(m+2,m+2), matriztoroide(m+2,m+2)
            !matrizes de pessoas que podem passar e matriz final para atualização
            integer :: passadores(m,m), matrizfim(m,m)
            !matriz de probabilidade de se pegar
            real :: prob(m,m)
            integer :: i, j, k, l


            call RANDOM_NUMBER(prob)
            passadores = 0
            matrizfim = matrix
            !não é mais um toroide
            matriztoroide(2:m+1,2:m+1) = matrix
            matrizt_toroide(2:m+1,2:m+1) = matriztimer
            
            !loop dos cálculos, vasculhamento das matrizes com influencias dos toroides
            !$OMP PARALLEL DO SIMD
            do j=2,m+1
                do i=2,m+1
                    !se a pessoa não tiver, tem uma checagem pra ver se tem gente perto que tenha
                    if ( matriztoroide(i,j) == 0 .and. sum(matriztoroide(i-1:i+1,j-1:j+1))  >= 1 ) then
                        !checagem se essas pessoas podem transmitir (timer) e somar os passadores 
                        do k = i-1,i+1
                            do l = j-1,j+1
                                if ( matrizt_toroide(k,l) >= 7 .and. matrizt_toroide(k,l) <= 25 ) then
                                    passadores(i-1,j-1) = passadores(i-1,j-1) + 1
                                end if
                            end do
                        end do
                        !se a probabilidade disse que a pessoa vai pegar, ela pega, tem o efeito
                        !de aglomeração também
                        if (passadores(i-1,j-1)*prob(i-1,j-1) >= 0.973 ) then
                            matrizfim(i-1,j-1) = 1
                        end if
                        !se a pessoa ta com sorte, ficou em casa ou seila, ela tem uma chance de não pegar
                        if ( prob(i-1,j-1) <= 0.25 ) then
                            matrizfim(i-1,j-1) = 2
                        end if
                    end if
                end do
            end do
            !$OMP END PARALLEL DO SIMD

        end function calculo_infeccao

        !função de atualização dos mortos
        function atualizacao_mortos(matriz,timer,mortos,m) result(mortes)
            USE OMP_LIB
            implicit none

            integer, intent(in) :: m
            integer :: i,j
            integer, intent(in) :: matriz(m,m)
            integer, intent(in) :: mortos(m,m), timer(m,m)
            integer :: mortes(m,m)
            real :: r(m,m)

            call RANDOM_NUMBER(r)

            mortes = mortos

            !se ela tiver, o timer apitar, e der a probabilidade, ela morre ou continua viva
            !$OMP PARALLEL DO SIMD
            do i = 1,m
                do j = 1,m
                    if ( matriz(j,i) == 1 .and. timer(j,i) >= 25 .and. mortos(j,i) == 0 ) then
                        if ( r(j,i) >= 0.99 ) then
                            mortes(j,i) = 1
                        else
                            mortes(j,i) = 2
                        end if
                    end if
                end do
            end do
            !$OMP END PARALLEL DO SIMD

        end function atualizacao_mortos

        !função de soma (1 para morto e 2 para vivo, porém com sequelas)
        function soma_mortos(matriz,m) result(somaa)
            USE OMP_LIB
            implicit none

            integer, intent(in) :: m
            integer, intent(in) :: matriz(m,m)
            integer :: somaa,i,j

            somaa = 0

            !cálculo da soma
            !$OMP PARALLEL DO SIMD REDUCTION(+:somaa)
            do i=1,m
                do j=1,m
                    if ( matriz(j,i) == 1 ) then
                        somaa = somaa+1
                    end if
                end do
            end do
            !$OMP END PARALLEL DO SIMD

        end function soma_mortos

        !soma dos infectados (1 para tem e 2 para deu sorte)
        function soma_infectados(matriz,m) result(somaa)
            USE OMP_LIB
            implicit none

            integer, intent(in) :: m
            integer, intent(in) :: matriz(m,m)
            integer :: somaa,i,j

            somaa = 0

            !soma propriamente dita
            !$OMP PARALLEL DO SIMD REDUCTION(+:somaa)
            do i=1,m
                do j=1,m
                    if ( matriz(j,i) == 1 ) then
                        somaa = somaa+1
                    end if
                end do
            end do
            !$OMP END PARALLEL DO SIMD

        end function soma_infectados

end module automatocelular