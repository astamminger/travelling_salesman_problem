module datatypes

  implicit none
  integer,         parameter    :: idp = selected_real_kind(15,307)
  !! Dieser Typ enthält die Information über die Koordinaten der 
  !! verwendeten Städte
  type city_t
    real(idp)                   :: rxcoordinate, rycoordinate
  end type city_t
  !! Enthält die Koodinateninformationen aller verwendeten Städte sowie 
  !! die Gesamtanzahl der Städte als Sammlung von einzelnen type city_t
  !! Datentypen
  type collection_t
    type(city_t),  allocatable  :: city(:)
    integer                     :: icity_count
  end type collection_t
  !! Enthält alle Informationen über den "Gen-Pool"
  type population_t
    integer                     :: ipop_count
    integer,       allocatable  :: ipop(:,:)
    integer,       allocatable  :: ipop_offspring(:,:)
    real(idp),     allocatable  :: rpop_fitness(:)
  end type population_t

  !! 
  type world_t
    type(population_t),  allocatable  :: pool(:)
  end type world_t

end module
