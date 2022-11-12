    
program HD
    use mpi

    real*8 :: a1,b1,a2,b2,a3,b3
    integer :: n1,n2,n3,i
    real*8,dimension(1000) :: xi,yi,zi,w1,w2,w3
    real*8, dimension(1000) :: local_wri, local_z
    real*8 :: GGS
    integer :: local_k
    integer :: my_rank
    integer :: ierr
    integer :: p
    real*8 :: total, summation




    call MPI_INIT(ierr)
 
        call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,  p, ierr)
            a1 = 0 
            a2 = 0
            a3 = 0
            b1 = 5
            b2 = 5
            b3 = 7
            n0 = 8
            call Get_data3(a1,b1,n0,my_rank,xi,w1)
            call Get_data4(a2,b2,n0,my_rank,yi,w2)
            call Get_data5(a3,b3,n0,my_rank,zi,w3)

            local_k = n0/p

            do i = 1,local_k
                local_z(i) = zi(i + my_rank*local_k)
                local_wri(i) = w3(i + my_rank*local_k)
                
            end do 

            summation = GGS(xi,yi,w1,w2,local_wri,local_z,local_k,n0,n0)


            call MPI_REDUCE( summation,  total, 1, MPI_DOUBLE ,MPI_SUM, 0, MPI_COMM_WORLD, ierr )

            if(my_rank .eq. 0) then
            print *, total 
            end if

    call MPI_FINALIZE(ierr)

end program 


    
    real*8 function GGS(x,y,w1,w2,w3,local_z,local_k,n1,n2)

    integer :: i,j,k,local_k
    integer :: n1,n2
    real*8 :: sum0, sum1 , sum2, integ
    real*8,dimension(1000) :: x,y,w1,w2,local_z,w3

    sum2 = 0.d0
    do k = 1, local_k 
        sum1 = 0.d0
        do i = 1, n2
            
            sum0 = 0.d0      
            do j = 1, n1
              
              integ = (x(j)*x(j)*x(j)+y(i)*y(i)*y(i)+local_z(k)*local_z(k)*local_z(k))
              sum0 = sum0 + integ*w1(j)

            enddo
            sum1 = sum1 + sum0*w2(i)
        enddo
        sum2 = sum2 + sum1*w3(k)
    enddo
    
    
    GGS = sum2*1
        
    return 
    end
      