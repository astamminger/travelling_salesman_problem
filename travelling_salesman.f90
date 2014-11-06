!! Programm zum Lösen des Travelling-Salesman Problem mit Hilfe eines 
!! genetischen Algorithmus. Nachkommen werden dabei mit Hilfe des
!! "Sequential Constructive Crossover" bzw. durch simple Mutation (Vertauschen
!! der Gene eines Individuums) erzeugt. Die Selektion ist durch eine
!! "Rank-based Roulette-Wheel Selection" realisiert.
!!
!! Info: Alle Subroutinen, die einen call auf die intrinsische Subroutine 
!! random_number() enthalten sind in eine !$CRITICAL Umgebung gesetzt, da
!! der parallele Aufruf von random_number() durch mehrere Threads anscheinend 
!! zu einer race condition führt, die das Programm u.U. extrem verlangsamt.
!
!
!
program travelling_salesman

  use datatypes

  implicit none

  character(len=2)         :: comment
  integer                  :: ierr, icity_num, i, j, k, l, thread_id, &
                              &min_loc(2), migs, min_pos
  integer,    allocatable  :: all_time_best(:), elite_genes(:,:,:), &
                              &mig_genes(:,:), best_ind(:,:,:)
  real(idp)                :: fit, tstart, tstop, rrndm_num, all_time_best_fit
  real(idp),  allocatable  :: costs(:,:), elite_fitness(:,:), mig_fitness(:), &
                              &min_fit(:), max_fit(:), mean_fit(:), &
                              &best_fit(:,:)
  type(collection_t)       :: col
  type(world_t)            :: world
  
  real(idp),  parameter    :: zero = 0.0d000
  real(idp),  parameter    :: prec = 1.0d-013
  real(idp),  parameter    :: pos_inf = huge(1.0d000)
 
  integer                  :: omp_get_thread_num
 
  !! Programm-Parameter
  !!
  !! Anzahl der verwendeten Individuen im "Gen-Pool"
  integer,    parameter  :: num_individuals = 250
  !! Anzahl an Crossover pro Generation
  integer,    parameter  :: n_cross = 10
  !! Mutationsrate (Für jede Position im Indiviuum)
  real(idp),  parameter  :: mutprob = 1.0d-002
  !! Anzahl an Generationen (Iterationsschritte)
  integer,    parameter  :: gen = 50000
  !! Anzahl paralleler Populationen (Anzahl verfügbarer Kerne)
  integer,    parameter  :: n_threads = 8
  !! Anzahl an Individuen die migriert werden
  integer,    parameter  :: n_mig = 1
  !! Anzahl an Generationen, nach denen Migration stattfindet
  integer,    parameter  :: n_gen_mig = 10000


  call cpu_time(tstart)
  call init_random_seed()

  call omp_set_num_threads(n_threads)
  call omp_set_dynamic(.false.)

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

  !! Mittlerer Fitnesswert vs. Iterationsschritt / Runtime
  open(unit=110,file='/tmp/mean_vs_it_step.dat',action='WRITE')
  !! Header
  write(110,'(A2,I23,A25)')comment, gen, 'lines'
  write(110,'(A2,A23,2A25)')comment, 'Iterationsschritt', 'Zeit',  '<Fitness>'

  !! Bester Fitnesswert vs. Iterationsschritt / Runtime
  open(unit=120,file='/tmp/best_vs_it_step.dat',action='WRITE')
  !! Header
  write(120,'(A2,I23,A25)')comment, gen, 'lines'
  write(120,'(A2,A23,2A25)')comment, 'Iterationsschritt', 'Zeit', &
                           &'Bester Fitnesswert'

  !! Initialisieren der Kostenmatrix
  allocate(costs(col%icity_count, col%icity_count))
  call init_cost_matrix(costs, col)


  allocate(elite_genes(n_threads, col%icity_count, num_individuals), &
          &elite_fitness(n_threads, num_individuals), &
          &mig_genes(col%icity_count, n_mig), &
          &mig_fitness(n_mig), min_fit(n_threads), max_fit(n_threads), &
          &mean_fit(n_threads), best_ind(n_threads, col%icity_count, 1), &
          &best_fit(n_threads, 1))


  allocate(world%pool(n_threads), all_time_best(col%icity_count))
  all_time_best = 0
  all_time_best_fit = pos_inf
  min_fit = pos_inf
  max_fit = zero
  mean_fit = zero

  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, thread_id) 

  thread_id = omp_get_thread_num() + 1

  allocate(world%pool(thread_id)%ipop(col%icity_count, num_individuals), &
          &world%pool(thread_id)%rpop_fitness(num_individuals), &
          &world%pool(thread_id)%ipop_offspring(col%icity_count, 2 * n_cross))
  world%pool(thread_id)%ipop = 0
  world%pool(thread_id)%ipop_offspring = 0
  world%pool(thread_id)%rpop_fitness = zero
  world%pool(thread_id)%ipop_count = num_individuals

  !! Erstellen der Startpopulation
  !$OMP CRITICAL
  call init_pool(world%pool(thread_id))
  !$OMP END CRITICAL

  do i = 1, gen, 1

    !! Sync Threads
    !$OMP BARRIER

    !! Speichern des besten Individuums
    if (i > 1) then
      min_pos = minloc(world%pool(thread_id)%rpop_fitness, 1)
      best_fit(thread_id, 1) = world%pool(thread_id)%rpop_fitness(min_pos)
      best_ind(thread_id, :, 1) = world%pool(thread_id)%ipop(:, min_pos)
    endif

    !$OMP CRITICAL
    call mutation(world%pool(thread_id), mutprob)
    call crossover(world%pool(thread_id), col, n_cross, mode='scr', &
                  &cost_matrix=costs, selection_pressure=real(1.10, idp)) 
    call rank_based_rw_selection(world%pool(thread_id), real(1.80,idp))
    call two_opt_optimization(world%pool(thread_id), col, costs)
    !$OMP END CRITICAL

    !! Nach n_gen_mig - Generationen schaffen es die insgesamt n_mig besten Gene
    !! in die anderen Pools zu migrieren
    if (mod(i, n_gen_mig) == 0) then

      elite_genes(thread_id, :, :) = 0
      elite_fitness(thread_id, :) = pos_inf

      elite_fitness(thread_id, :) = world%pool(thread_id)%rpop_fitness(:)
      elite_genes(thread_id, :, :) = world%pool(thread_id)%ipop(:, :)
      !! Warten bis alle Threads ihre Daten geschrieben haben
      !$OMP BARRIER

      !! Die n_mig besten Kandidaten auswählen. Jedes Gen nur einmal, keine 
      !! Duplikate, da sonst bei einem bereits stark konvergierten Pool dann 
      !! immer nur Gene von diesem Pool migriert werden würden auch wenn der
      !! Pool nicht gegen das globale Optimum konvergiert ist.
      !$OMP MASTER
      k = 1
      mig_fitness = zero
      mig_genes = 0
      do
        min_loc = 0
        min_loc = minloc(elite_fitness(:, :))
        if (k == 1) then 
          mig_fitness(k) = minval(elite_fitness)
          mig_genes(:, k) = elite_genes(min_loc(1), :, min_loc(2))
          k = k + 1
        else
          !! Ist der neue minimale Wert nicht gleich dem vorherigen Wert wird
          !! er akzeptiert
          if ( abs(minval(elite_fitness) - mig_fitness(k - 1)) > prec ) then
            mig_fitness(k) = minval(elite_fitness)
            mig_genes(:, k) = elite_genes(min_loc(1), :, min_loc(2))
            k = k + 1
          endif
        endif
        elite_fitness(min_loc(1), min_loc(2)) = pos_inf
        !! Sind n_mig unterschiedliche Gene gefunden oder nicht genügend
        !! unterschiedliche Gene vorhanden:  Verlassen der Schleife.
        if (k > n_mig) exit
        if (abs(minval(elite_fitness) - pos_inf) <= prec) exit
      enddo  
      !! Werden keine n_mig unterschiedliche Gene gefunden, werden die 
      !! verbleibenden Stellen mit dem besten gefundenen Gen besetzt.
      if (minval(mig_genes(:,:)) == 0) then
        do l = 1, n_mig, 1
          if (minval(mig_genes(:, l)) == 0) then
            mig_genes(:, l) = mig_genes(:, 1)
            mig_fitness(l) = mig_fitness(1)
          endif
        enddo
      endif
      !$OMP END MASTER

      !! Synchronisieren von Master und den anderen Threads, damit alle Daten
      !! vollständig vorhanden sind
      !$OMP BARRIER

      !! Ersetzen der n_mig schlechtesten Gene in jedem Pool durch 
      !! die ausgewählten besten Gene aus allen Pools
      migs = num_individuals - n_mig + 1
      call order_ascending(world%pool(thread_id))
      world%pool(thread_id)%rpop_fitness(migs:) = mig_fitness(:)
      world%pool(thread_id)%ipop(:,migs:) = mig_genes(:,:)
    endif

    !! Ist das beste Individuum schlechter als das beste Individuum der 
    !! Vorgängergeneration wird es durch den Vorgänger ersetzt.
    if (i > 1) then
      min_pos = minloc(world%pool(thread_id)%rpop_fitness, 1)
      if (minval(world%pool(thread_id)%rpop_fitness) &
         &> best_fit(thread_id, 1)) then
        world%pool(thread_id)%rpop_fitness(min_pos) = best_fit(thread_id, 1)
        world%pool(thread_id)%ipop(:, min_pos) = best_ind(thread_id, :, 1)
      endif
    endif

    !! Ausgeben von Informationen über den momentanten Status
    if (mod(i, 10) == 0 .or. i == 1) then
      min_fit(thread_id) = minval(world%pool(thread_id)%rpop_fitness)
      max_fit(thread_id) = maxval(world%pool(thread_id)%rpop_fitness)
      mean_fit(thread_id) = mean_fitness(world%pool(thread_id))
      !$OMP BARRIER
      !$OMP MASTER
      write(*,*)' '
      write(*,'(A12,I10)')'Generation: ', i
      write(*,'(A8,3A25)')'TaskID: ', '<Fitness>', 'Min(Fitness)', 'Max(Fitness)'
      do j = 1, n_threads, 1
        write(*,'(I8,3ES25.15)')j - 1, mean_fit(j), min_fit(j), max_fit(j)
      enddo
      !! Ausgabe in File: Gesamtbester Fitnesswert und Gesamtmittelwert
      call cpu_time(tstop)
      write(110, '(I25,2ES25.15)') i, (tstop - tstart) / n_threads, &
                                  &sum(mean_fit) / n_threads
      write(120, '(I25,2ES25.15)') i, (tstop - tstart) / n_threads, &
                                  &minval(min_fit)
      !! Ausgabe der Gesamtbesten Route
      call dump_best_to_file(world%pool(minloc(min_fit, 1)), col, &
                            &filename='/tmp/best_tour.dat')
      !$OMP END MASTER
    endif

  enddo

  !$OMP END PARALLEL

  close(110)
  close(120)

  !! Freigeben des verwendeten Speichers
  call apocalypse()


contains


  !! Wählt zufälliges Individuum aus dem Pool aus und entfernt Kreuzungen aus
  !! der Tour bzw. fügt Kreuzungen in die Tour ein sofern dies die Fitness
  !! verbessert.
  subroutine two_opt_optimization(population, collection, cost_matrix)
    type(population_t),       intent(inout)  :: population
    type(collection_t),       intent(inout)  :: collection
    real(idp),  allocatable,  intent(inout)  :: cost_matrix(:,:)
    real(idp),  allocatable                  :: mod_fit(:)
    real(idp)                                :: rrndm_num, change_fit
    integer,    allocatable                  :: path(:), mod_path(:,:)
    integer                                  :: i, j, n_city, ind

    n_city = size(population%ipop(:, 1))
    allocate(path(n_city), mod_path(n_city, 1), mod_fit(1))
    mod_path = 0
    path = 0
    change_fit = zero

    !! Wähle zufälliges Individuum
    call random_number(rrndm_num)
    ind = ceiling(rrndm_num * size(population%ipop(1,:)))
    path = population%ipop(:, ind)
    fit = population%rpop_fitness(ind)
    mod_path(:,1) = path(:)

    do i = 1, n_city - 2, 1
      do j = i + 2, n_city - 1, 1
        !! Änderung der Fitness beim Tausch der Städte ermitteln  
        change_fit = cost_matrix(path(i),path(j)) + &
                     &cost_matrix(path(i+1),path(j+1)) - &
                     &cost_matrix(path(i),path(i+1)) - &
                     &cost_matrix(path(j),path(j+1))
        if (change_fit < zero) then
          mod_path(j,1) = path(i + 1)
          mod_path(i + 1,1) = path(j)
          path(:) = mod_path(:,1)
        endif  
      enddo
    enddo
  
    mod_path(:, 1) = path(:)
    call get_path_length(mod_path, mod_fit, collection) 

    !! Überschreiben des schlechtesten Individuums sofern ein besserer Wert
    !! gefunden wurde
    if (fit /= mod_fit(1)) then
      population%ipop(:, maxloc(population%rpop_fitness, 1)) = mod_path(:,1)
      population%rpop_fitness(maxloc(population%rpop_fitness, 1)) = mod_fit(1)
    endif

    deallocate(path, mod_path, mod_fit)
  end subroutine

  !! Initialisieren der Ausgangspopulation. Dazu wird das Array pop_array 
  !! zunächst geordnet mit Zahlen von 1 bis pop_sizex gefüllt und anschließend 
  !! die Reihenfolge der Zahlen pop_sizex-mal zuällig vertauscht. 
  subroutine init_pool(population)
    type(population_t)  :: population
    real(idp)           :: rrndm_num
    integer             :: i, n_city, n_pop

    n_city = size(population%ipop(:,1))
    n_pop = size(population%ipop(1,:))

    !! Geordnetes Füllen des Arrays
    do i = 1, n_city, 1
      population%ipop(i,:) = i
    enddo
  
    !! Reihenfolge der Zahlen zufällig vertauschen (willkürlich n_city-mal) 
    do i = 1, n_pop, 1
      call swap_digits(population%ipop(:,i), n_city)
    enddo
  end subroutine

  !! Subroutine zum Mutieren eines Individuums. Jede Position der Tour hat
  !! dabei die Wahrscheinlichkeit mut_probability mit einer anderen zufälligen
  !! Position der Tour vertauscht zu werden. 
  subroutine mutation(population, mut_probability)
    type(population_t),  intent(inout)  :: population 
    real(idp),           intent(in)     :: mut_probability
    real(idp)                           :: rrndm_num
    integer                             :: i, j, n_pop, pos1, pos2, temp, &
                                           &n_city

    n_pop = size(population%ipop(1,:))
    n_city = size(population%ipop(:,1))

    do i = 1, n_pop, 1
      do j = 1, n_city, 1
        rrndm_num = zero
        pos1 = 0
        pos2 = 0

        !! Wenn rrndm_num =< mut_probability wird die Mutation ausgeführt
        call random_number(rrndm_num)
        if (rrndm_num <= mut_probability) then
          call random_number(rrndm_num)
          pos1 = j
          pos2 = ceiling(rrndm_num * n_city)

          !! Vertauschen der beiden Positionen
          temp = population%ipop(pos1, i)
          population%ipop(pos1, i) = population%ipop(pos2, i)
          population%ipop(pos2, i) = temp
          temp = 0
        endif
      enddo
    enddo
  end subroutine

  !! Crossover Subroutine. Es werden zunächst die beiden Eltern über eine 
  !! rangbasierte Roulette-Rad Selektion gewählt. Im Anschluss daran wird
  !! die eigentliche Subroutine für den Crossover anhand des in mode gegebenen
  !! Werts aufgerufen und der Crossover ausgeführt. (mode = standard oder scr)
  subroutine crossover(population, collection, crossover_count, mode, &
                      &cost_matrix, selection_pressure)
    type(population_t),       intent(inout)  :: population
    type(collection_t),       intent(inout)  :: collection
    real(idp),  allocatable,  optional :: cost_matrix(:,:)
    real(idp),  allocatable            :: temp_fit(:), sorted_fit(:), prob(:)
    real(idp)                          :: rrndm_num, norm, selection_pressure
    integer,    allocatable            :: parents(:,:), childs(:,:), &
                                          &sorted(:,:), temp(:,:)
    integer                            :: i, n_city, n_pop, crossover_count, &
                                          &maxfit_pos, j, k
    character(len=*)                   :: mode

    n_city = size(population%ipop(:,1))
    n_pop = size(population%ipop(1,:))

    allocate(parents(n_city, 2), childs(n_city, 2))
    parents = 0
    childs = 0

    !! Prüfe ob genug Individuuen für die gewünschte Anzahl an Crossovers 
    !! vorhanden sind (Mindestens n_pop >= 2*crossover_count)
    if (n_pop - 2 * crossover_count < 0) then
      write(*,*)'Number of desired crossovers exceeding number of individuals!'
      stop
    elseif (n_pop - 2 * crossover_count == 1) then
        write(*,*)'Reduce number of crossovers per step by one!'
    endif

    !! Zuweisen der temporären Arrays
    allocate(temp(n_city, n_pop), temp_fit(n_pop), sorted(n_city, n_pop), &
            &sorted_fit(n_pop), prob(n_pop))

    temp = 0
    sorted = 0
    temp_fit = zero
    sorted_fit = zero
    prob = zero
    norm = zero

    !! Übergeben der Population an temp(:,:)
    temp = population%ipop

    call get_path_length(temp, temp_fit, col) 

    !! Sortieren der Einträge in aufsteigender Reigenfolge. Damit entspricht
    !! der Rang des Individuums gerade seiner Position in sorted(:,:)
    do i = 1, size(temp_fit), 1
      maxfit_pos = maxloc(temp_fit, 1)
      sorted(:, i) = temp(:, maxfit_pos)
      sorted_fit(i) = temp_fit(maxfit_pos)
      temp_fit(maxfit_pos) = real(-1, idp)
    enddo


    !! Berechnung der Wahrscheinlichkeitsanteile der verschiedenen Ränge
    !! für die anschließende Roulette Auswahl. Je höher der Rang, desto
    !! größer ist auch der Segmentanteil im "Roulette-Rad".
    do i = 1, size(prob), 1
      prob(i) = real(2, idp) - selection_pressure + real(2, idp) * (& 
                &selection_pressure - real(1, idp)) * (real(i, idp) &
                & - real(1, idp)) / (real(n_pop, idp)  - real(1, idp))
      norm = norm + prob(i)

      !! Bereiche = Segmentanteil am "Roulette-Rad" der einzelnen Ränge
      !! bestimmen
      if (i > 1) then
        prob(i) = prob(i) + prob(i - 1)
      endif
    enddo

    prob = prob / norm

    do i = 1, crossover_count, 1

      !! Würfle 2 Eltern anhand der rbrw-Selection 
      do j = 1, 2, 1
        call random_number(rrndm_num)
         do k = 1, size(prob), 1
           if (prob(k) > rrndm_num) then
             parents(:,j) = sorted(:, k)
             exit
           endif
         enddo   
      enddo
  
      if (mode == 'standard') then
        call standard_crossover(parents, childs)
      elseif (mode == 'scr') then
        call sequential_constructive_crossover(collection, cost_matrix, &
                                               &parents, childs(:, 1)) 
      else
        write(*,'(A)')'Mode not supported!'
        stop
      endif

      !! Ersetzen der Eltern mit den neu erzeugten Kindern
      if (mode == 'standard') then
        population%ipop_offspring(:, 2 * i - 1) = childs(:, 1)
        population%ipop_offspring(:, 2 * i) = childs(:, 2)
      elseif (mode == 'scr') then
        population%ipop_offspring(:, i) = childs(:, 1)
      endif
    enddo

    deallocate(temp, temp_fit, sorted, sorted_fit, prob, parents, childs)
  end subroutine

  !! Standard Crossover wählt zufällig zwei Position und überträgt die Sequenz
  !! zwischen diesen beiden Positionen von Elternteil 1 auf Kind 2 und von
  !! Elternteil 2 auf Kind 1. Die verbleibenden Positionen von Kind 1 werden
  !! dann mit den Genen des Elternteil 1 befüllt die nicht bereits in der
  !! übertragenen Sequenz enthalten waren (Kind 2 analog)
  subroutine standard_crossover(parents, childs)
    integer    ::  parents(:,:), childs(:,:)
    integer    :: j, k, l, n_city, ipos(2)
    real(idp)  :: rrndm_num
    logical    :: duplicate

    n_city = size(parents(:, 1))

    !! Erzeugen der Kinder aus den gewählten Eltern
    do
      call random_number(rrndm_num)
      ipos(1) = ceiling(n_city * rrndm_num)
      call random_number(rrndm_num)
      ipos(2) = ceiling(n_city * rrndm_num)
      if (ipos(1) /= ipos(2)) exit
    enddo
    !! Sortieren der Position, so dass ipos(1) der kleinere Wert ist
    if (ipos(2) < ipos(1)) then
      k = ipos(1)
      ipos(1) = ipos(2)
      ipos(2) = k
    endif

    !! Ausgewählten Bereich der Elternteile auf die Kinder übertragen
    do j = ipos(1), ipos(2), 1
      childs(j, 2) = parents(j, 1)
      childs(j, 1) = parents(j, 2)
    enddo
    !! Verbleibende Plätze mit Städten aus dem jeweils anderen Elternteil 
    !! füllen. Entspricht eine Zahl des Elternteils einer Zahl des bereits
    !! geschriebenen Bereichs wird diese übergangen.
    do l = 1, 2, 1
      do j = 1, n_city, 1
        do k = 1, n_city, 1
          duplicate = .false.
          duplicate = check_duplicate(childs(:, l), parents(k, l))
          if (j < ipos(1) .or. j > ipos(2)) then
            if (duplicate .eqv. .false.) then
              childs(j, l) = parents(k, l)
              exit
            endif
          endif
        enddo
      enddo
    enddo
  end subroutine

  !! Sequential constructive crossover erzeugt ein Kind aus zwei Eltern und 
  !! gewichtet dazu die Wahl der Knoten mit einer Kostenmatrix, die die Kosten
  !! für jede direkte Verbindung von zwei Knotenpunken enthält. Aus den ersten 
  !! und zweiten Elternteil wird je eine Verbindung zwischen zwei Knoten erzeugt
  !! und im Anschluss die mit den geringsten Kosten gewählt und auf das Kind
  !! übertragen
  subroutine sequential_constructive_crossover(collection, rcost_matrix, &
                                               &parents, child)
    type(collection_t)       :: collection
    real(idp)                :: rcost_matrix(:,:)
    integer                  :: parents(:,:), child(:)
    integer                  :: node(2), node_pos(2)
    real(idp)                :: cost_node(2)
    integer                  :: i, j, n_city

    n_city = collection%icity_count 

    node = 0
    node_pos = 0
    child = 0
    do i = 1, n_city, 1
      !! Erster Knoten des Kindes entspricht erstem Knoten des 
      !! ersten Elternteils
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
              &.eqv. .false.) then
            node(1) = parents(node_pos(1) + 1, 1)
          !! Knoten bereits vorhanden, dann wähle den nächsten Knoten aus
          !! parents(:, 1) der noch nicht in child(:) enthalten ist.
          else
            do j = 1, n_city, 1
              if (check_duplicate(child(:), parents(j, 1)) .eqv. .false.) then
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
            if (check_duplicate(child(:), parents(j, 1)) .eqv. .false.) then
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
              &.eqv. .false.) then
            node(2) = parents(node_pos(2) + 1, 2)
          !! Knoten bereits vorhanden, dann wähle den nächsten Knoten aus
          !! parents(:, 1) der noch nicht in childs(:, 1) enthalten ist.
          else
            do j = 1, n_city, 1
              if (check_duplicate(child(:), parents(j, 2)) .eqv. .false.) then
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
            if (check_duplicate(child(:), parents(j, 2)) .eqv. .false.) then
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
  end subroutine

  !! Initialisieren der Kosten-Matrix. Die Matrix enthält die Kosten (Weglänge)
  !! für alle direkten Verbindungen zwischen den Knoten. Die Kosten entsprechen
  !! hier einfach der 2D euklidischen Norm des Vektors von Stadt i zu Stadt j.
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
 
  !! Berechnen der Wegstrecke (Fitness) für ein Individuum des Pools. Die
  !! berechnete Weglänge wird dann entsprechend der Position in pop_array(:,i)
  !! an die zugehörige Position i im  Array fitness(:) geschrieben.
  subroutine get_path_length(population, fitness, city_collection)
    integer,    allocatable,  intent(inout)  :: population(:, :)
    real(idp),  allocatable,  intent(inout)  :: fitness(:)
    type(collection_t),        intent(in)    :: city_collection
    integer                                  :: i, j, n_pop, n_city, ipos, &
                                                &iposp1
    real(idp)                                :: rxn, rxnp1, ryn, rynp1     

    n_city = city_collection%icity_count
    n_pop = size(population(1, :))
   
    rxn = zero
    rxnp1 = zero
    ryn = zero
    rynp1 = zero
    fitness = zero

    !! Berechnung der Weglänge (Fitness) des i-ten Individuums
    do i = 1, n_pop, 1
      do j = 1, n_city - 1, 1
        ipos = population(j,i)
        iposp1 = population(j + 1, i)
        rxn = city_collection%city(ipos)%rxcoordinate
        rxnp1 = city_collection%city(iposp1)%rxcoordinate
        ryn = city_collection%city(ipos)%rycoordinate
        rynp1 = city_collection%city(iposp1)%rycoordinate

        fitness(i) = fitness(i) + sqrt((rxn - rxnp1) * (rxn - rxnp1) + &
                     &(ryn - rynp1) * (ryn - rynp1))
      enddo

      !! Zurück zum Anfangspunkt
      ipos = population(n_city,i)
      iposp1 = population(1, i)
      rxn = city_collection%city(ipos)%rxcoordinate
      rxnp1 = city_collection%city(iposp1)%rxcoordinate
      ryn = city_collection%city(ipos)%rycoordinate
      rynp1 = city_collection%city(iposp1)%rycoordinate
      fitness(i) = fitness(i) + sqrt((rxn - rxnp1) * (rxn - rxnp1) + &
                   &(ryn - rynp1) * (ryn - rynp1))
    enddo
  end subroutine

  !! Wählt Individuen aus modified_population gemäß ihrer Fitness aus und 
  !! schreibt die Individuen mit der besten Fitness in population. Der Wert
  !! von survival_prob gibt dabei die Überlebenswahrscheinlichkeit von 
  !! Individuen mit schlechterer Fitness an. Nach der Auswahl der besten
  !! Gene werden die schlechtesten Gene durch die vorhandenen noch schlechteren
  !! Gene ersetzt, sofern eine gewürfelte Zufallszahl kleiner ist als der 
  !! angegebene Wert von survival_prob,
  subroutine selection(population, survival_prob)
    type(population_t),  intent(inout)  :: population
    integer,       allocatable          :: temp(:, :), sorted(:, :)
    real(idp),     allocatable          :: temp_fit(:), sorted_fit(:)
    real(idp)                           :: rrndm_num, survival_prob
    integer                             :: i, minfit_pos, n_pop, n_city, n_off,&
                                           & N, j

    n_pop = population%ipop_count
    n_city = size(population%ipop(:, 1))
    !! Abzählen der tatsächlich vorhandenen Kinder
    do i = 1, size(population%ipop_offspring(1, :)), 1
      if (population%ipop_offspring(1, i) == 0) then
        n_off = i - 1
        exit
      elseif(i == size(population%ipop_offspring(1, :))) then
        n_off = i
      endif
    enddo
    N = n_pop + n_off

    !! Allocate temporäres Array mit der Doppelten Anzahl an Individuen
    !! falls alle Individuen in modified_population verändert wurden
    allocate(temp(n_city, N), temp_fit(N), sorted(n_city, N), sorted_fit(N))

    temp = 0
    sorted = 0
    temp_fit = zero
    sorted_fit = zero

    !! Zusammenfassen von gesamter Population und Offspring in einem
    !! Array
    temp(:, 1:n_pop) = population%ipop(:, 1:n_pop)
    temp(:, (n_pop + 1):) = population%ipop_offspring(:, 1:n_off)

    call get_path_length(temp, temp_fit, col) 
    
    !! Sortieren der Einträge in absteigender Reigenfolge. D.h bester zuerst
    !! schlechtester am Ende.
    do i = 1, size(temp_fit), 1
      minfit_pos = minloc(temp_fit, 1)
      sorted(:, i) = temp(:, minfit_pos)
      sorted_fit(i) = temp_fit(minfit_pos)
      temp_fit(minfit_pos) = pos_inf
    enddo

    do i = 1, n_pop, 1
      population%ipop(:, i) = sorted(:, i)
      population%rpop_fitness(i) = sorted_fit(i)
    enddo

    !! Ersetze die schlechtesten Gene mit der Wahrscheinlichkeit
    !! survival_prob durch die noch schlechteren Gene im Pool
    do i = 1, n_off, 1
      call random_number(rrndm_num)
      if (rrndm_num < survival_prob) then
        population%ipop(:, n_pop - n_off + i) = sorted(:, n_pop + i)
        population%rpop_fitness(n_pop - n_off + i) = sorted_fit(n_pop + i)
      endif 
    enddo

    deallocate(temp, temp_fit, sorted, sorted_fit)
  end subroutine

  !! Verbesserte Selektion, die die Wahrscheinlichkeit für das
  !! Überleben eines Individuums nicht an den Fitnesswert koppelt sondern an
  !! den Rang (entspricht der Position im Array) den das Individuum anhand 
  !! seiner Fitness zugewiesen bekommt. In diesem Fall gilt je niedriger der 
  !! Rang desto schlechter das Indiviuun (hoher Fitnesswert).
  !! Anhand des Rangs wird dann über eine lineare Funktion die
  !! Wahrscheinlichkeit bestimmt mit der das Individuum im Roulette-Rad
  !! vertreten ist. (Diese Wahrscheinlichkeit kann über den Selektionsdruck
  !! selection_pressure kontrolliert werden. selection_pressure entspricht 
  !! dabei gerade der Wahrscheinlichkeit mit der das fitteste Individuum 
  !! im Vergleich zu allen anderen Individuen gewählt wird. Die 
  !! Wahrscheinlichkeit, das das schlechteste Indiviuum gewählt wird 
  !! entspricht gerade selecion_pressure - 1) 
  subroutine rank_based_rw_selection(population, selection_pressure)
    type(population_t),  intent(inout)  :: population
    integer,       allocatable          :: temp(:, :), sorted(:, :)
    real(idp),     allocatable          :: temp_fit(:), sorted_fit(:), prob(:)
    real(idp),           intent(in)     :: selection_pressure
    real(idp)                           :: rrndm_num, norm
    integer                             :: i, maxfit_pos, n_pop, n_city, n_off,&
                                           & N, j

    n_pop = population%ipop_count
    n_city = size(population%ipop(:, 1))
    !! Abzählen der tatsächlich vorhandenen Kinder
    do i = 1, size(population%ipop_offspring(1, :)), 1
      if (population%ipop_offspring(1, i) == 0) then
        n_off = i - 1
        exit
      elseif(i == size(population%ipop_offspring(1, :))) then
        n_off = size(population%ipop_offspring(1, :))
      endif
    enddo
    N = n_pop + n_off

    !! Allocate temporäres Array mit der Doppelten Anzahl an Individuen
    !! falls alle Individuen in modified_population verändert wurden
    allocate(temp(n_city, N), temp_fit(N), sorted(n_city, N), sorted_fit(N), &
            &prob(N))

    temp = 0
    sorted = 0
    temp_fit = zero
    sorted_fit = zero
    prob = zero
    norm = zero

    !! Zusammenfassen von gesamter Population und Offspring in einem
    !! Array
    temp(:, 1:n_pop) = population%ipop(:, 1:n_pop)
    temp(:, (n_pop + 1):) = population%ipop_offspring(:, 1:n_off)


    call get_path_length(temp, temp_fit, col) 
    
    !! Sortieren der Einträge in aufsteigender Reigenfolge. Damit entspricht
    !! der Rang des Individuums gerade seiner Position in sorted(:,:)
    do i = 1, size(temp_fit), 1
      maxfit_pos = maxloc(temp_fit, 1)
      sorted(:, i) = temp(:, maxfit_pos)
      sorted_fit(i) = temp_fit(maxfit_pos)
      temp_fit(maxfit_pos) = real(-1, idp)
    enddo


    !! Berechnung der Wahrscheinlichkeitsanteile der verschiedenen Ränge
    !! für die anschließende Roulette Auswahl. Je höher der Rang, desto
    !! größer ist auch der Segmentanteil im "Roulette-Rad". 
    do i = 1, size(prob), 1
      prob(i) = real(2, idp) - selection_pressure + real(2, idp) * (& 
                &selection_pressure - real(1, idp)) * (real(i, idp) &
                & - real(1, idp)) / (real(N, idp)  - real(1, idp))
      norm = norm + prob(i)

      !! Bereiche = Segmentanteil am "Roulette-Rad" der einzelnen Ränge
      !! bestimmen
      if (i > 1) then
        prob(i) = prob(i) + prob(i - 1)
      endif
    enddo

    prob = prob / norm

    !! Würfle n_pop-mal. Das Individuum in dessen Segment die Zufallszahl liegt
    !! wird dann in die neue Population geschrieben
    do j = 1, n_pop, 1
      call random_number(rrndm_num)
       do i = 1, N, 1
         if (prob(i) > rrndm_num) then
           population%ipop(:, j) = sorted(:, i)
           population%rpop_fitness(j) = sorted_fit(i)
           exit
         endif
       enddo
     enddo

    deallocate(temp, temp_fit, sorted, sorted_fit, prob)
  end subroutine

  !! Subroutine zum Sortieren der Population in aufsteigender Reihenfolge. 
  !! Das heißt auf Platz 1 befindet sich dann das Gen mit der niedrigsten 
  !! Fitness und auf dem letzten Platz das Gen mit der höchsten Fitness.
  !! (Hier: niedrigste Fitness = beste Fitness da Algorithmus als 
  !! Minimierungsproblem geschrieben wurde)
  subroutine order_ascending(population)
    type(population_t),  intent(inout)  :: population
    integer,  allocatable               :: temp(:, :), sorted(:, :)
    real(idp),  allocatable             :: temp_fit(:), sorted_fit(:)
    integer                             :: i, N, n_city, minfit_pos 

    N = size(population%ipop(1, :))
    n_city = size(population%ipop(:, 1))

    !! Allocate temporäres Array mit der Doppelten Anzahl an Individuen
    !! falls alle Individuen in modified_population verändert wurden
    allocate(temp(n_city, N), temp_fit(N), sorted(n_city, N), sorted_fit(N))

    temp = 0
    sorted = 0
    temp_fit = zero
    sorted_fit = zero

    temp = population%ipop
    temp_fit = population%rpop_fitness
    
    !! Sortieren der Einträge in aufsteigender Reigenfolge. Damit entspricht
    !! der Rang des Individuums gerade seiner Position in sorted(:,:)
    do i = 1, size(temp_fit), 1
      minfit_pos = minloc(temp_fit, 1)
      sorted(:, i) = temp(:, minfit_pos)
      sorted_fit(i) = temp_fit(minfit_pos)
      temp_fit(minfit_pos) = pos_inf
    enddo

    population%ipop = sorted
    population%rpop_fitness = sorted_fit

    deallocate(temp, temp_fit, sorted, sorted_fit)
  end subroutine

  !! Schreibt die derzeit beste Route in population in eine in filename
  !! spezifizierte Datei. Es wird sowohl die Route als auch die zu den
  !! Städten zugehörigen Koordinaten geschrieben.
  subroutine dump_best_to_file(population, city_collection, filename)
    type(population_t),  intent(in)  :: population
    type(collection_t),  intent(in)  :: city_collection
    character(len=*),     intent(in)  :: filename
    integer                          :: i, ipos, icity

    open(unit=100, file=trim(adjustl(filename)), action='WRITE')
    write(100,'(A2,A23,2A25)')'# ', 'Stadt', 'x-Koordinate', 'y-Koordinate'
    ipos = minloc(population%rpop_fitness, 1)

    do i = 1, city_collection%icity_count, 1
      icity = population%ipop(i, ipos)
      write(100,'(I25,2ES25.15)')icity, &
                                 &city_collection%city(icity)%rxcoordinate, &
                                 &city_collection%city(icity)%rycoordinate
    enddo
    !! Schreibe Startpunkt an (icity_count + 1)-te Stelle im File
    !! (Erleichtert plotten der Tour mit gnuplot).
    icity = population%ipop(1, ipos)
    write(100,'(I25,2ES25.15)')icity, &
                               &city_collection%city(icity)%rxcoordinate, &
                               &city_collection%city(icity)%rycoordinate

    close(100)
  end subroutine

  !! Subroutine zum freigeben des gesamten, verwendeten Speichers
  subroutine apocalypse()

    deallocate(col%city)
    deallocate(world%pool)
    deallocate(costs)
    deallocate(all_time_best)
    deallocate(elite_genes)
    deallocate(elite_fitness)
    deallocate(mig_genes)
    deallocate(mig_fitness)
    deallocate(best_fit)
    deallocate(best_ind)
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
  
  !! Funktion um zu prüfen ob eine bestimmte Stadt int_value in array(:) 
  !! bereits enthalten ist.
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

  !! Funktion zum ermittelt der Position einer bestimmten Stadt (node) 
  !! in array(:).
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

  !! Funkion zur Berechnung des Mittelwertes der Fitness der in population
  !! gespeicherten Population.
  function mean_fitness(population)
    real(idp)                       :: mean_fitness
    type(population_t),  intent(in) :: population
    integer                         :: i
    real(idp)                       :: fitness_sum

    fitness_sum = zero
    do i = 1, size(population%rpop_fitness), 1
      fitness_sum = fitness_sum + population%rpop_fitness(i)
    enddo

    mean_fitness = fitness_sum / size(population%rpop_fitness)
  end function

  !! Initialisiert einen Seed für den Zufallszahlengenerator anhand der 
  !! momentanen Systemzeit.
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
