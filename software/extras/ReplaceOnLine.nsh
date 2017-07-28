!macro ReplaceOnLine SEARCH_TEXT REPLACEMENT LINE SOURCE_FILE
  Push "${SEARCH_TEXT}"
  Push "${REPLACEMENT}"
  Push "${LINE}"
  Push "${SOURCE_FILE}"
  Call ROnL
!macroend

Function ROnL
Exch $R0 ;file
Exch
Exch $R1 ;line
Exch 2
Exch $R2 ;replace with
Exch 2
Exch 3
Exch $R3 ;string to replace
Exch 3
Push $R4
Push $R5
Push $R6
Push $R7
Push $R8
Push $R9
Push $9
 FileOpen $R4 $R0 r
 GetTempFileName $R5
 FileOpen $R6 $R5 w
Top:
  FileRead $R4 $R7
  IntOp $R8 $R8 + 1
  StrCmp $R8 $R1 +3
  FileWrite $R6 $R7
  Goto Top
   StrLen $9 $R3
Loop_Top:
   StrCpy $R8 0
Loop:
   IntOp $R8 $R8 - 1
   StrCpy $R9 $R7 $9 $R8
   StrCmp $R9 "" Finish
   StrCmp $R9 $R3 0 Loop
    StrCpy $R9 $R7 $R8
    IntOp $R8 $R8 + $9
    StrCpy $R7 $R7 "" $R8
    StrCpy $R7 $R9$R2$R7
detailprint $R7
    FileWrite $R6 $R7
    Goto Loop_Top
Finish:
  ClearErrors
  FileRead $R4 $R7
  IfErrors +3
  FileWrite $R6 $R7
  Goto Finish
FileClose $R4
FileClose $R6
SetDetailsPrint none
Delete $R0
Rename $R5 $R0
SetDetailsPrint both
Pop $9
Pop $R9
Pop $R8
Pop $R7
Pop $R6
Pop $R5
Pop $R4
Pop $R3
Pop $R2
Pop $R1
Pop $R0
FunctionEnd