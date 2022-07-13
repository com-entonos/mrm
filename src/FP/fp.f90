module FP_m
  use, intrinsic :: ISO_C_Binding, only: C_int, C_ptr, C_NULL_ptr, C_DOUBLE, C_ASSOCIATED

  implicit none
  private
  public :: new_FP, advance_FP, delete_FP, copy_FP, null_FP, FP_s

  type FP_s
    private
    type(C_ptr) :: object = C_NULL_ptr
  end type
  interface
    function C_FP_new () result(this) bind(C,name="CFP_new")
      import
      type(C_ptr) :: this
    end function
    function C_FP_copy (this) result(that) bind(C,name="CFP_copy")
      import
      type(C_ptr),value :: this
      type(C_ptr) :: that
    end function
    subroutine C_FP_advance (this,mx,my,mz,kB_T_over_mu,kB_Tc_over_mu,dt) bind(C,name="CFP_advance")
      import
      real(C_DOUBLE)::mx,my,mz,kB_T_over_mu,kB_Tc_over_mu,dt
      type(C_ptr),value::this
    end subroutine
    subroutine C_FP_delete (this) bind(C,name="CFP_delete")
      import
      type(C_ptr),value::this
    end subroutine

  end interface
  interface new_FP
    module procedure CFP_new
  end interface new_FP
  interface copy_FP
    module procedure CFP_copy
  end interface copy_FP
  interface advance_FP
    module procedure CFP_advance
  end interface advance_FP
  interface delete_FP
    module procedure CFP_delete
  end interface delete_FP
contains
! Fortran wrapper routines to interface C wrappers
  subroutine CFP_new(this)
    type(FP_s), intent(out) :: this
    this%object = C_FP_new()
  end subroutine
  subroutine CFP_delete(this)
    type(FP_s), intent(inout) :: this
    call C_FP_delete(this%object)
    this%object = C_NULL_ptr
  end subroutine
  subroutine CFP_copy(this, that)
    type(FP_s), intent(in) :: this
    type(FP_s), intent(out) :: that
    that%object = C_FP_copy(this%object)
  end subroutine
  subroutine CFP_advance(this,mx,my,mz,kB_T_over_mu,kB_Tc_over_mu,dt)
    type(FP_s), intent(in) :: this
    real(C_DOUBLE),intent(inout)::mx,my,mz
    real(C_DOUBLE),intent(in)::kB_T_over_mu,kB_Tc_over_mu,dt
    call C_FP_advance(this%object,mx,my,mz,kB_T_over_mu,kB_Tc_over_mu,dt)
  end subroutine
  function null_FP(this) result (null_p)
    type(FP_s), intent(in) :: this
    logical::null_p
    null_p=.not.c_associated(this%object)
  end function

end module
