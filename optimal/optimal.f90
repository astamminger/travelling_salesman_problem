program optimal
 
  use datatypes
 
  implicit none

  type(collection_t)       :: col
  integer                  :: i, N, ierr, icity_num
  character(len=2)         :: comment
  integer,  allocatable    :: combinations(:,:), work(:)
  real(idp),  allocatable  :: path_length(:)
  integer                  :: k, i1, i2, i3, i4, i5, i6, i7, i8, i9, i10

  open(unit=100, file='./cities.dat', action='read')
  read(100,*,iostat=ierr) comment, col%icity_count
  if (ierr == 0) then
    allocate(col%city(col%icity_count))
  else
    write(*,*)'Error'
    stop
  endif
  !! Koordinaten einlesen
  i = 1
  do
    read(100,*,iostat=ierr) icity_num, col%city(i)%rxcoordinate, &
                           &col%city(i)%rycoordinate
    if (ierr == 0) then
      i = i + 1
    endif
    if (ierr < 0 ) exit           !! ierr < 0 -> Dateiende
    if (i > col%icity_count) exit 
  enddo
  close(100)

icity_num = 9
  allocate(work(icity_num), path_length(faculty(icity_num)), &
          &combinations(icity_num, faculty(icity_num)))
  path_length = 0
  work = 0
  N = icity_num
  k = 1
  do i1 = 1, N, 1
  do i2 = 1, N, 1
  do i3 = 1, N, 1
  do i4 = 1, N, 1
  do i5 = 1, N, 1
  do i6 = 1, N, 1
  do i7 = 1, N, 1
  do i8 = 1, N, 1
  do i9 = 1, N, 1
  work(1) = i1
  work(2) = i2
  work(3) = i3
  work(4) = i4
  work(5) = i5
  work(6) = i6
  work(7) = i7
  work(8) = i8
  work(9) = i9

  if (unique(work) == .true.) then
    combinations(:,k) = work(:)
    path_length(k) = pathl(work, col)
  if (mod(k, 10000) == 0) then
write(*,'(I7,A,I7)')k,' von ', faculty(icity_num)
  elseif (k == faculty(icity_num)) then
write(*,'(I7,A,I7)')k, ' von ', k
  endif
    k = k + 1
  endif
  enddo
  enddo
  enddo
  enddo
  enddo
  enddo
  enddo
  enddo
  enddo
write(*,*)' '
write(*,'(A,ES25.15)')'Minimale Pfadlänge: ', minval(path_length)
write(*,'(A,9I3)')'Minimaler Pfad: ', combinations(:, minloc(path_length, 1))

contains

  function faculty(n)
    integer               :: faculty
    integer,  intent(in)  :: n
    integer               :: i

    faculty = 1
    do i = 2, n, 1
      faculty = faculty * i 
    enddo
  end function

  function unique(work)
    logical               :: unique
    integer,  intent(in)  :: work(:) 
    integer               :: i, j, abort

    unique = .true.
    abort = 0
    do i = 1, size(work), 1
      do j = i + 1, size(work), 1
        if (work(i) == work(j)) then
          unique = .false.
          abort = 1
          exit
        endif  
      enddo
      if (abort == 1) then
        exit
      endif
    enddo

  end function

  function pathl(combination, collection)
    real(idp) :: pathl
    integer,  intent(in) :: combination(:)
    type(collection_t), intent(in) :: collection
    real(idp)   :: x, xp1, y, yp1
    integer     :: i, cit_no, cit_ne
     

    pathl = 0
    do i = 1, size(combination) - 1, 1
      cit_no = combination(i)
      cit_ne = combination(i + 1)

      x = collection%city(cit_no)%rxcoordinate
      y = collection%city(cit_no)%rycoordinate
      xp1 = collection%city(cit_ne)%rxcoordinate
      yp1 = collection%city(cit_ne)%rycoordinate

      pathl = pathl + sqrt((x - xp1) * (x - xp1) + (y - yp1) * (y - yp1))
    enddo

  end function
   
       



end program
