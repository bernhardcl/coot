Function TextCompare
	!define TextCompare `!insertmacro TextCompareCall`
 
	!macro TextCompareCall _FILE1 _FILE2 _OPTION _FUNC
		Push $0
		Push `${_FILE1}`
		Push `${_FILE2}`
		Push `${_OPTION}`
		GetFunctionAddress $0 `${_FUNC}`
		Push `$0`
		Call TextCompare
		Pop $0
	!macroend
 
	Exch $3
	Exch
	Exch $2
	Exch
	Exch 2
	Exch $1
	Exch 2
	Exch 3
	Exch $0
	Exch 3
	Push $4
	Push $5
	Push $6
	Push $7
	Push $8
	Push $9
	ClearErrors
 
	IfFileExists $0 0 error
	IfFileExists $1 0 error
	StrCmp $2 'FastDiff' +5
	StrCmp $2 'FastEqual' +4
	StrCmp $2 'SlowDiff' +3
	StrCmp $2 'SlowEqual' +2
	goto error
 
	FileOpen $4 $0 r
	IfErrors error
	FileOpen $5 $1 r
	IfErrors error
	SetDetailsPrint textonly
 
	StrCpy $6 0
	StrCpy $8 0
 
	nextline:
	StrCmp $4 '' fast
	IntOp $8 $8 + 1
	FileRead $4 $9
	IfErrors 0 +4
	FileClose $4
	StrCpy $4 ''
	StrCmp $5 '' end
	StrCmp $2 'FastDiff' fast
	StrCmp $2 'FastEqual' fast slow
 
	fast:
	StrCmp $5 '' call
	IntOp $6 $6 + 1
	FileRead $5 $7
	IfErrors 0 +5
	FileClose $5
	StrCpy $5 ''
	StrCmp $4 '' end
	StrCmp $2 'FastDiff' call close
	StrCmp $2 'FastDiff' 0 +2
	StrCmp $7 $9 nextline call
	StrCmp $7 $9 call nextline
 
	slow:
	StrCmp $4 '' close
	StrCpy $6 ''
	DetailPrint '$8. $9'
	FileSeek $5 0
 
	slownext:
	FileRead $5 $7
	IfErrors 0 +2
	StrCmp $2 'SlowDiff' call nextline
	StrCmp $2 'SlowDiff' 0 +2
	StrCmp $7 $9 nextline slownext
	IntOp $6 $6 + 1
	StrCmp $7 $9 0 slownext
 
	call:
	Push $2
	Push $3
	Push $4
	Push $5
	Push $6
	Push $7
	Push $8
	Push $9
	Call $3
	Pop $0
	Pop $9
	Pop $8
	Pop $7
	Pop $6
	Pop $5
	Pop $4
	Pop $3
	Pop $2
	StrCmp $0 'StopTextCompare' 0 nextline
 
	close:
	FileClose $4
	FileClose $5
	goto end
 
	error:
	SetErrors
 
	end:
	SetDetailsPrint both
	Pop $9
	Pop $8
	Pop $7
	Pop $6
	Pop $5
	Pop $4
	Pop $3
	Pop $2
	Pop $1
	Pop $0
FunctionEnd
