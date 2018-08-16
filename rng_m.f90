module rng_m
  use precision_m
  use MersenneTwister, only: ignpoi

  interface rand_num
     module procedure rand_num_double
     module procedure rand_num_array
  end interface

  public :: rand_num,seedrnd

contains

  ! -------------------------------------------------------

  subroutine seedrnd()
    ! DESCRIPTION:
    !	Seeds the random number generator.
    Use MersenneTwister, only: sgrnd
    Implicit None
    integer :: seed, clock

    !        Integer :: sd(1)

    call system_clock(count=clock)

    seed = clock + 37 

    call sgrnd(seed)

    !        sd(1) = seed
    !        Call Random_Seed(PUT=sd)
  end subroutine seedrnd

  ! -------------------------------------------------------

  subroutine rand_num_double(rnd)
    Use MersenneTwister, only: grnd
    implicit none
    real(wp), intent(out) :: rnd
    
    rnd = grnd()
    !        Call Random_Number(Rnd)
  end subroutine rand_num_double

  subroutine rand_num_array(rnd)
    Use MersenneTwister, only: grnd
    implicit none
    real(wp), intent(out) :: rnd(:)
    integer :: num,i

    num = size(rnd)
    do i=1,num
       rnd(i) = grnd()
    end do
  end subroutine rand_num_array


  ! -------------------------------------------------------

  real(wp) function IRnd(top)
    ! DESCRIPTION:
    !	Returns a uniform random integer in [1,top].
    Use MersenneTwister, only: igrnd
    Implicit None

    Integer, Intent(IN) :: top
    !        Real :: num

    IRnd = IAnd(igrnd(), top-1) + 1

    !        Call Random_Number(num)
    !        IRnd = Int(num * Real(top) - 0.5E0) + 1
  end function IRnd

end module rng_m
