program gencity

  implicit none

  integer,    parameter    :: idp = selected_real_kind(15,307)
  integer                  :: i
  character(len=12)         :: path_to_file = './cities.dat'
  integer,    parameter    :: icity_count = 50   !! Anzahl Städte
  real(idp),  parameter    :: rrange = 1d001    !! Bereich der Koordinaten
                                                !! x,y in [0,rrange]
  real(idp)                :: rrandom_number    !! Zufallszahl
  real(idp),  allocatable  :: rcoordinates(:,:) !! Koordinaten 

  !! Initialisieren des zufälligen Seeds
  call init_random_seed()

  !! Erzeugen der Zufallskoordinaten
  allocate(rcoordinates(2,icity_count))

  open(unit=100, file=path_to_file, action='write')
write(100,'(A2,I3,A8)')'# ', icity_count, ' Städte'
write(100,'(A2,A5,2A25)')'# ', 'Stadt', 'x-Koordinate', 'y-Koordinate'

  do i = 1, icity_count

    !! x-Koordinate
!    call random_number(rrandom_number) 
!    rcoordinates(1,i) = ceiling(rrandom_number * rrange)
   
    !! y-Koordinate
!    call random_number(rrandom_number)
!    rcoordinates(2,i) = ceiling(rrandom_number * rrange)


!write(100,'(I6,2ES25.15)')i, real(50,idp)*(real(1,idp)+cos(real(2*3.14,idp)/real(icity_count,idp)*real(i-1,idp))), real(50,idp)*(real(1,idp)+sin(real(2*3.14,idp)/real(icity_count,idp)*real(i-1,idp)))


  enddo
 

!    write(100,'(A2,I3,A8)')'# ', icity_count, ' Städte'
!    write(100,'(A2,A5,2A25)')'# ', 'Stadt', 'x-Koordinate', 'y-Koordinate'

!    do i = 1, icity_count, 1
!      write(100,'(I6,2ES25.15)')i, rcoordinates(1,i), rcoordinates(2,i)
!    enddo

  close(100)
 
  deallocate(rcoordinates)


contains


  subroutine init_random_seed()
    integer :: i, n, clock
    integer, dimension(:), allocatable :: seed

    call random_seed(size = n)
    allocate(seed(n))

    call system_clock(count=clock)

    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    call random_seed(put = seed)

    deallocate(seed)
  end subroutine


end program
