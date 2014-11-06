program modified_sequential_constructive_crossover

  use datatypes

  implicit none
  
  integer,    allocatable  :: parents(:,:), child(:)
  integer                  :: node(2), node_pos(2)
  real(idp),  allocatable  :: rcost_matrix(:,:)
  real(idp)                :: cost_node(2)
  integer                  :: ierr, i, j, k, icity_num
  character(len=2)         :: comment = '# '
  type(collection_t)       :: col
  
  integer                  :: n_city
  real(idp),  parameter    :: pos_inf = huge(1.0d000)
  real(idp),  parameter    :: zero = 0.0d000 


  !! Einlesen der Städe und ihrer Koordinaten. Dabei wird zunächst die Anzahl
  !! der verwendeten Städte ermittelt und anschließend die Koordinaten
  !! in die entsprechenden Arrays geschrieben
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
allocate(parents(col%icity_count, 2), child(col%icity_count))
n_city = col%icity_count
parents = 0
child = 0
!! Nur hier zum testen !!!!!!!!!!!!!!!!!!!  
  do i = 1, n_city, 1                    !
    parents(i, 1) = i                    !
    parents(i, 2) = i                    !
  enddo                                  !
!! Nur hier zum testen !!!!!!!!!!!!!!!!!!!
  call swap_digits(parents(:,1), 100)    !
  call swap_digits(parents(:,2), 100)    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call sequential_constructive_crossover(col, parents, child)

contains

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

!! Sequential constructive crossover erzeugt ein Kind aus zwei Eltern und 
!! gewichtet dazu die Wahl der Knoten mit einer Kostenmatrix, die die Kosten
!! für jede direkte Verbindung von zwei Knotenpunken enthält. Aus den ersten 
!! und zweiten Elternteil wird je eine Verbindung zwischen zwei Knoten erzeugt
!! und im Anschluss die mit den geringsten Kosten gewählt und auf das Kind
!! übertragen
!!
!! collection: type(collection_t)
!! parents: 2D Array parents(n_city, 2) enthält die beiden Elternteile
!! child: 1D Array child(n_city) das später das erzeugte Kind enthält
subroutine sequential_constructive_crossover(collection, parents, child)
  type(collection_t)       :: collection
  integer                  :: parents(:,:), child(:)
  integer                  :: node(2), node_pos(2)
  real(idp),  allocatable  :: rcost_matrix(:,:)
  real(idp)                :: cost_node(2)
  integer                  :: i, j, n_city

  n_city = collection%icity_count

  !! Initialisieren der Kostenmatrix, die die Kosten für jede direkte 
  !! Verbindung zwischen zwei Knoten enthält
  allocate(rcost_matrix(n_city, n_city))
  rcost_matrix = zero
 
  call init_cost_matrix(rcost_matrix, collection) 


  node = 0
  node_pos = 0
  do i = 1, n_city, 1
    !! Erster Knoten des Kindes entspricht erstem Knoten des ersten Elternteils
    if (i == 1) then
      node(1) = parents(i,1)
      child(i) = node(1)
    
    !! Ermitteln der Nachfolgenden Knoten
    else

      !! Bestimmten der Position des ersten Knotens in parents(:,1) und 
      !! parents(:,2)
      node_pos(1) = node_position(parents(:, 1), node(1))
      node_pos(2) = node_position(parents(:, 2), node(1))

      !! Bestimmen von node(1) aus parents(:, 1)
      !! Liegt der Knoten nicht am Ende des Strings, prüfen ob Nachfolgeknoten
      !! Bereits in child enthalten ist. 
      if (node_pos(1) /= n_city) then
        if (check_duplicate(child(:), parents(node_pos(1) + 1, 1)) &
            &== .false.) then
          node(1) = parents(node_pos(1) + 1, 1)
        !! Knoten bereits vorhanden, dann wähle den nächsten Knoten aus
        !! parents(:, 1) der noch nicht in child(:) enthalten ist.
        else
          do j = 1, n_city, 1
            if (check_duplicate(child(:), parents(j, 1)) == .false.) then
              node_pos(1) = j
              node(1) = parents(j, 1)
              exit
            endif
          enddo
        endif
      !! Liegt der Knoten am Ende des Strings, dann wähle den nächsten Knoten
      !! aus parents(:, 1) der noch nicht in child(:) enthalten ist.
      elseif (node_pos(1) == n_city) then 
        do j = 1, n_city, 1
          if (check_duplicate(child(:), parents(j, 1)) == .false.) then
            node_pos(1) = j
            node(1) = parents(j, 1)
            exit
          endif
        enddo
      endif

      !! Bestimmen von node(2) aus parents(:, 2) analog zur Bestimmung
      !! von node(1). 
      if (node_pos(2) /= n_city) then
        if (check_duplicate(child(:), parents(node_pos(2) + 1, 2)) &
            &== .false.) then
          node(2) = parents(node_pos(2) + 1, 2)
        !! Knoten bereits vorhanden, dann wähle den nächsten Knoten aus
        !! parents(:, 1) der noch nicht in childs(:, 1) enthalten ist.
        else
          do j = 1, n_city, 1
            if (check_duplicate(child(:), parents(j, 2)) == .false.) then
              node_pos(2) = j
              node(2) = parents(j, 2)
              exit
            endif
          enddo
        endif
      !! Liegt der Knoten am Ende des Strings, dann wähle den nächsten Knoten
      !! aus parents(:, 1) der noch nicht in childs(:, 1) enthalten ist.
      elseif (node_pos(2) == n_city) then 
        do j = 1, n_city, 1
          if (check_duplicate(child(:), parents(j, 2)) == .false.) then
            node_pos(2) = j
            node(2) = parents(j, 2)
            exit
          endif
        enddo
      endif
  
      !! Jetzt den besten Nachfolge-Knoten anhand der Kosten-Matrix bestimmen.
      !! Sind beide Knoten gleich oder besitzen die gleichen Kosten wird 
      !! willkürlich Node1 als Nachfolger festgelegt. Sonst der mit den 
      !! niedrigsten Kosten 
      if (node(1) /= node(2)) then
        cost_node(1) = rcost_matrix(child(i), node(1))
        cost_node(2) = rcost_matrix(child(i), node(2))
        if (cost_node(1) > cost_node(2)) then
          child(i) = node(2)
          node(1) = node(2)
        elseif (cost_node(2) > cost_node(1)) then
          child(i) = node(1)
        elseif (cost_node(1) == cost_node(2)) then
          child(i) = node(1)
        endif
        elseif (node(1) == node(2)) then
          child(i) = node(1)
      endif  
    endif
  enddo
write(*,'(20I3)')parents(:,1)  
write(*,'(ES25.15)')pathl(parents(:,1), collection)
write(*,'(20I3)')parents(:,2)  
write(*,'(ES25.15)')pathl(parents(:,2), collection)
write(*,'(20I3)')child(:)
write(*,'(ES25.15)')pathl(child , collection)

  
  deallocate(rcost_matrix)
end subroutine

  subroutine init_cost_matrix(cost_matrix, collection)
    real(idp),           intent(inout)  :: cost_matrix(:,:)
    type(collection_t),  intent(in)     :: collection
    integer                             :: i, j
    real(idp)                           :: rxi, rxj, ryi, ryj

    do i = 1, size(cost_matrix(1,:)), 1
      do j = 1, size(cost_matrix(1,:)), 1
        rxi = collection%city(i)%rxcoordinate
        rxj = collection%city(j)%rxcoordinate
        ryi = collection%city(i)%rycoordinate
        ryj = collection%city(j)%rycoordinate
        
        if (j /= i) then    
          cost_matrix(i, j) = sqrt((rxi - rxj) * (rxi - rxj) &
                                  &+(ryi - ryj) * (ryi - ryj))
        !! Punkt a nach Punkt a ist keine sinnvolle Lösung. Da diese Lösungen 
        !! nicht berücksichtigt werden sollen werden die Kosten auf "unendlich"
        !! gesetzt
        else
          cost_matrix(i, j) = pos_inf
        endif
      enddo
    enddo
   
  end subroutine 
 
  !! Wählt n_swap mal zwei zufällige Werte aus array(:) aus und vertauscht
  !! deren Plätze im Array. 
  subroutine swap_digits(array, n_swap)
    integer,  intent(inout)  :: array(:)
    integer,  intent(in)     :: n_swap
    integer                  :: i, ipos(2), iwork(2), iarr_size
    real(idp)                :: rrndm_num
    iarr_size = size(array(:))

    i = 1
    do
      call random_number(rrndm_num)
      ipos(1) = ceiling(iarr_size * rrndm_num)
      call random_number(rrndm_num)
      ipos(2) = ceiling(iarr_size * rrndm_num)

      if (i > n_swap) then
        exit
      endif
      if (ipos(1) /= ipos(2)) then
        iwork(1) = array(ipos(1))
        iwork(2) = array(ipos(2))
     
        !! Tauschen der Werte
        array(ipos(1)) = iwork(2)
        array(ipos(2)) = iwork(1)

        i = i + 1
      endif  
    enddo
    
  end subroutine
 
  function node_position(array, node)

    integer                  :: node_position
    integer,  intent(inout)  :: array(:)
    integer,  intent(in)     :: node
    integer                  :: i 
     
    do i = 1, size(array), 1
      if (array(i) == node) then
        node_position = i
      endif
    enddo

  end function

  !! Prüft ob int_value in array(:) bereits enthalten ist
  function check_duplicate(array, int_value)

    logical               :: check_duplicate 
    integer,  intent(in)  :: array(:)
    integer,  intent(in)  :: int_value
    integer               :: i
   
    check_duplicate = .false. 
    do i = 1, size(array), 1
      if (array(i) == int_value) then
        check_duplicate = .true.
        exit
      endif
    enddo
   
  end function

end program
