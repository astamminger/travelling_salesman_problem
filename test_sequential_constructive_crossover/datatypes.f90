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
  !! Da bei gegebener Überlebensrate Mutation zu starken Fluktuationen
  !! führen können, auch wenn das Ergebnis bereits konvergiert ist,
  !! wird bei Mutationen keine Überlebenswahrscheinlichkeit verwendet.
  !! Um die Mutierten Individuen in der Selektion zu filtern wird daher 
  !! bei Mutation mut_flags(i) = .true. gesetzt.
  type population_t
    integer                     :: ipop_count
    integer,       allocatable  :: ipop(:,:)
    real(idp),     allocatable  :: rpop_fitness(:)
    logical,       allocatable  :: mut_flags(:)
  end type population_t

end module
