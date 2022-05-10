!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This module contains subroutines and interfaces to deal with arrays

MODULE Dealing_with_arrays
implicit none

! this interface finds by itself which of the two subroutine to use depending on the dimensions of the array passed:
interface Find_in_array ! search cheaking one by one
    module procedure Find_in_1D_array
    module procedure Find_in_2D_array
end interface Find_in_array
! this interface finds by itself which of the two subroutine to use depending on the dimensions of the array passed:
interface Find_in_array_monoton ! search with bisection method
    module procedure Find_in_monotonous_1D_array
    module procedure Find_in_monotonous_2D_array
end interface Find_in_array_monoton
! this interface finds by itself which of the two subroutine to use depending on the dimensions of the array passed:
interface Linear_approx
    module procedure Linear_approx_1d
    module procedure Linear_approx_2d
end interface Linear_approx

!private  ! hides items not listed on public statement 
public :: Find_in_array, Find_in_array_monoton, Linear_approx


contains

subroutine Linear_approx_1d(Array_x, Array_y, In_val, Value1, El1, El2)
   REAL(8), dimension(:), INTENT(in) :: Array_x ! grid for the array
   REAL(8), dimension(:), INTENT(in) :: Array_y ! in this array make an approximation
   real(8), INTENT(in) :: In_val    ! this is the value
   !integer, intent(in) :: Number    ! element of the array
   real(8), intent(in), optional :: El1, El2    ! where to start approximation from, if needed
   real(8), intent(out) :: Value1   ! output, approximated value
   real(8) el_one
   integer Number
   
   call Find_in_array_monoton(Array_x, In_val, Number)    ! find the closest value in the array to a given one
   
   if (Number .EQ. 1) then
       if (present(El1)) then  ! the starting points are known, use them:
          if (present(El2)) then ! the starting points are known, use them:
             Value1 = El2+(Array_y(Number)-El2)/(Array_x(Number)-El1)*(In_val - El1)
          else ! only the X-starting point is known, Y is not, use some assumptions:
             if (size(Array_y) .GE. 2) then ! array is long enough to assume something:
                if (Array_y(2) .GT. Array_y(1)) then ! it is locally decreasing array, assume starting point "infinity"
                    Value1 = 1d21
                else    ! it is locally increasing, assume starting point 'zero'
                    Value1 = Array_y(Number)/(Array_x(Number)-El1)*(In_val - El1)
                endif
             else   ! array is too short, no assumption can be made, just make it equal to the first value:
                Value1 = El1
             endif
          endif
       else ! no starting points are present, nothing to assume, use first point as approximation:
          Value1 = Array_y(1) ! [A] total mean free path
       endif
   else
      if (Array_y(Number-1) .GT. 1d20) then ! if it starts from infinity, approximate as 'infinity'
         Value1 = Array_y(Number-1)
      else  ! if it's normal array, just interpolate:
         Value1 = Array_y(Number-1)+(Array_y(Number)-Array_y(Number-1))/(Array_x(Number)-Array_x(Number-1))*(In_val - Array_x(Number-1))
      endif
   endif
end subroutine Linear_approx_1d

subroutine Linear_approx_2d(Array, In_val, Value1, El1, El2)
   REAL(8), dimension(:,:), INTENT(in) :: Array ! in this array make an approximation
   real(8), INTENT(in) :: In_val    ! this is the value
   !integer, intent(in) :: Number    ! element of the array
   real(8), intent(in), optional :: El1, El2    ! where to start approximation from, if needed
   real(8), intent(out) :: Value1   ! output, approximated value
   real(8) el_one
   integer Number
   
   call Find_in_array_monoton(Array, In_val, 1, Number)    ! find the closest value in the array to a given one
   
   if (Number .EQ. 1) then
       if (present(El1)) then  ! the starting points are known, use them:
          if (present(El2)) then ! the starting points are known, use them:
             Value1 = El2+(Array(2,Number)-El2)/(Array(1,Number)-El1)*(In_val - El1)
          else ! only the X-starting point is known, Y is not, use some assumptions:
             if (size(Array,2) .GE. 2) then ! array is long enough to assume something:
                if (Array(2,1) .GT. Array(2,2)) then ! it is locally decreasing array, assume starting point "infinity"
                    Value1 = 1d21
                else    ! it is locally increasing, assume starting point 'zero'
                    Value1 = Array(2,Number)/(Array(1,Number)-El1)*(In_val - El1)
                endif
             else   ! array is too short, no assumption can be made, just make it equal to the first value:
                Value1 = El1
             endif
          endif
       else ! no starting points are present, nothing to assume, use first point as approximation:
          Value1 = Array(2,1) ! [A] total mean free path
       endif
   else
      if (Array(2,Number-1) .GT. 1d20) then ! if it starts from infinity, approximate as 'infinity'
         Value1 = Array(2,Number-1)
      else  ! if it's normal array, just interpolate:
         Value1 = Array(2,Number-1)+(Array(2,Number)-Array(2,Number-1))/(Array(1,Number)-Array(1,Number-1))*(In_val - Array(1,Number-1))
      endif
   endif
end subroutine Linear_approx_2d

subroutine Find_in_1D_array(Array, Value, Number)
   REAL(8), dimension(:), INTENT(in) :: Array ! in which we are looking for the Value
   REAL(8), INTENT(in) :: Value   ! to be found in the array as near as possible
   integer, INTENT(out) :: Number ! number of the element which we are looking for 
   integer i
   i = 1
   do while (Array(i) .LT. Value)
      i = i + 1
   enddo
   Number = i
end subroutine Find_in_1D_array

subroutine Find_in_2D_array(Array, Value, Indx, Number)
   REAL(8), dimension(:,:), INTENT(in) :: Array ! in which we are looking for the Value
   REAL(8), INTENT(in) :: Value   ! to be found in the array as near as possible
   integer, INTENT(in) :: Indx    ! index of the array, showing in which colonm we search
   integer, INTENT(out) :: Number ! number of the element which we are looking for 
   integer i
   i = 1
   do while (Array(Indx,i) .LT. Value)
      i = i + 1
   enddo
   Number = i
end subroutine Find_in_2D_array

subroutine Find_in_monotonous_1D_array(Array, Value0, Number)
   REAL(8), dimension(:), INTENT(in) :: Array ! in which we are looking for the Value
   REAL(8), INTENT(in) :: Value0   ! to be found in the array as near as possible
   integer, INTENT(out) :: Number ! number of the element which we are looking for 
   integer i, N, i_cur, i_1, i_2, coun
   real(8) temp_val, val_1, val_2
   
   N = size(Array)
   i_1 = 1
   val_1 = Array(i_1)
   i_2 = N
   val_2 = Array(i_2)
   i_cur = FLOOR((i_1+i_2)/2.0)
   temp_val = Array(i_cur)
   if (isnan(Value0)) then
        print*, 'The subroutine Find_in_monotonous_1D_array'
        print*, 'cannot proceed, the value of Value0 is', Value0
        write(*, '(f,f,f,f)') Value0, Array(i_cur), Array(i_1), Array(i_2)
        pause 'STOPPED WORKING...'
   else
       if (Value0 .LT. Array(1)) then ! it's the first value, no need to search
           i_cur = 0
       else if (Value0 .GE. Array(N)) then ! it's the last value, no need to search
           i_cur = N-1
       else
           coun = 0
           do ! until the Value is in between Array(i_cur) and Array(i_cur+1) => we found i_cur
                if ((Value0 .GE. Array(i_cur)) .AND. (Value0 .LE. Array(i_cur+1))) exit ! when the Value is in between Array(i_cur) and Array(i_cur+1) => we found i_cur
                if (temp_val .LE. Value0) then
                   i_1 = i_cur
                   val_1 = Array(i_1)
                   i_cur = FLOOR((i_1+i_2)/2.0)
                   temp_val = Array(i_cur)
                else
                   i_2 = i_cur
                   !val_2 = Array(i_2)
                   val_2 = temp_val
                   i_cur = FLOOR((i_1+i_2)/2.0)
                   temp_val = Array(i_cur)
                endif
                coun = coun + 1
                if (coun .GT. 1e3) then
                    print*, 'PROBLEM WITH CONVERGANCE IN'
                    print*, 'Find_in_monotonous_1D_array', coun
                    write(*, '(f,f,f,f)') Value0, Array(i_cur), Array(i_1), Array(i_2)
                    pause 'STOPPED WORKING...'
                endif
           enddo
       endif
   endif    ! isnan
   Number = i_cur+1
end subroutine Find_in_monotonous_1D_array

subroutine Find_in_monotonous_2D_array(Array, Value0, Indx, Number)
   REAL(8), dimension(:,:), INTENT(in) :: Array ! in which we are looking for the Value
   REAL(8), INTENT(in) :: Value0   ! to be found in the array as near as possible
   integer, INTENT(in) :: Indx    ! index of the array, showing in which colonm we search
   integer, INTENT(out) :: Number ! number of the element which we are looking for 
   integer i, N, i_cur, i_1, i_2, coun
   real(8) temp_val, val_1, val_2
   
   N = size(Array,2)
   i_1 = 1
   val_1 = Array(Indx,i_1)
   i_2 = N
   val_2 = Array(Indx,i_2)
   i_cur = FLOOR((i_1+i_2)/2.0)
   temp_val = Array(Indx,i_cur)
   
   if (isnan(Value0)) then
        print*, 'The subroutine Find_in_monotonous_2D_array'
        print*, 'cannot proceed, the value of Value0 is', Value0
        write(*, '(f,f,f,f)') Value0, Array(Indx,i_cur), Array(Indx,i_1), Array(Indx,i_2)
        pause 'STOPPED WORKING...'
   else
       if (Value0 .LT. Array(Indx,1)) then ! it's the first value, no need to search
           i_cur = 0
       else if (Value0 .GE. Array(Indx,N)) then ! it's the last value, no need to search
           i_cur = N-1
       else
           coun = 0
           do ! until the Value is in between Array(i_cur) and Array(i_cur+1) => we found i_cur
                if ((Value0 .GE. Array(Indx,i_cur)) .AND. (Value0 .LE. Array(Indx,i_cur+1))) exit ! when the Value is in between Array(i_cur) and Array(i_cur+1) => we found i_cur
                if (temp_val .LE. Value0) then
                   i_1 = i_cur
                   val_1 = temp_val
                   i_cur = FLOOR((i_1+i_2)/2.0)
                   temp_val = Array(Indx,i_cur)
                else
                   i_2 = i_cur
                   val_2 = temp_val
                   i_cur = FLOOR((i_1+i_2)/2.0)
                   temp_val = Array(Indx,i_cur)
                endif
                
                coun = coun + 1
                if (coun .GT. 1e3) then
                    print*, 'PROBLEM WITH CONVERGANCE IN'
                    print*, 'Find_in_monotonous_2D_array', coun
                    write(*, '(f,f,f,f)') Value0, Array(Indx,i_cur), Array(Indx,i_1), Array(Indx,i_2)
                    pause 'STOPPED WORKING...'
                endif
           enddo
       endif
   endif    ! isnan
   Number = i_cur+1
end subroutine Find_in_monotonous_2D_array


END MODULE Dealing_with_arrays