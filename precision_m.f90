module precision_m
  implicit none

  integer,parameter :: wp=selected_real_kind(15,307)
  private
  public :: wp
end module precision_m
