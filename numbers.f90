!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Module: numbers
!Purpose: contains declarations for parameters relating to 
!         variable declarations. Primarily used for dealing
!         with the KIND atrribute.
!Author: Reese Horton, 19456155
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module numbers

  !kind parameter for real variables we want (15 decimal places, 10^+-307 exponent, old fortran double)
  !Used in all other modules
  !Syntax for assigning values is x_dp e.g. 2.0_dp gives 2.0 to the precision dp
  integer, parameter:: dp = selected_real_kind(15,307)  



end module numbers



