program cmbspec
  use healpix_types
  use params
  use time_mod
  implicit none

  ! Initialize time grids
  call initialize_time_mod
  write(*,*) "Hello world!"

  ! Output to file desired quantities here

end program cmbspec
