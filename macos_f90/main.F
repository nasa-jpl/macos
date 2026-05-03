	Function GetNumToken(Val)
        Implicit none

        Character(len=220) :: Val
        Integer :: GetNumToken, nbc,c
        Logical :: inTok

        Do nbc=1,220
          if (Val(nbc:nbc) /= ' ') exit
        End Do
        If (nbc .GE. 220) Then
          GetNumToken=0
          return
        End If

        inTok=.false.; GetNumToken=0
        Do c=nbc,220
          if (.not. inTok .and. Val(c:c) /= ' ') then
            GetNumToken = GetNumToken + 1; inTok=.true.
          else if (inTok .and. Val(c:c)==' ') then
            inTok=.false.
          end if
        End Do
	print*,' ** GetNumToken =',GetNumToken
        End Function GetNumToken

	Program main
	Character(len=220) :: msg
	Integer :: nTok,GetNumToken

	msg = "333 5.6d0 78d-3 hello, how are you?"
	nTok = GetNumToken(msg)
	print*,' ** nTok =',nTok
	End Program main
