PROGRAM MAIN
    IMPLICIT NONE

CONTAINS

SUBROUTINE test
    IMPLICIT NONE
    type linkedl(na)
        integer(4), len :: na
        real, dimension(na) :: list
        type(linkedl(na)), pointer :: next
    end type linkedl

    type(linkedl(100)), pointer :: head=>null()
    type(linkedl(100)), pointer :: tail=>null()
    type(linkedl(100)), pointer :: tmp=>null()

    print *, associated(head)
    print *, associated(tail)
    print *, associated(tmp)
END SUBROUTINE
END PROGRAM
