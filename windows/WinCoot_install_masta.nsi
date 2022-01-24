; Script generated by the HM NIS Edit Script Wizard.

; a temporary fixed hack to not make guile, will be variable at some point
; Could make with a guile version
!define WITH_GUILE
; Do not make a guile version
!undef WITH_GUILE

;for bat change
!include "ReplaceOnLine.nsh"
!include "AdvReplaceInFile.nsh"
!include "TextCompare.nsh"

; for logicals
!include "LogicLib.nsh"

; for command line options
!include FileFunc.nsh
!insertmacro GetParameters
!insertmacro GetOptions

; for removal of installer
!include "StrStr.nsh"

; to detect windows version
!include "WinVer.nsh"

!ifndef src_dir
!define src_dir "C:\msys64\home\bernhard\autobuild\MINGW64_NT-6.1-bernie-pre-release-gtk2-shared"
!endif
; pre setting of Coot version
!include "${src_dir}\coot-version"

; HM NIS Edit Wizard helper defines
!define PRODUCT_NAME "WinCoot"
!define PRODUCT_VERSION "${WinCootVersion}"
!define PRODUCT_PUBLISHER "Bernhard Lohkamp & Paul Emsley"
!define PRODUCT_WEB_SITE "http://bernhardcl.github.io/coot"
!define PRODUCT_WEB_SITE_2 "http://www2.mrc-lmb.cam.ac.uk/personal/pemsley/coot/"
!define PRODUCT_DIR_REGKEY "Software\Microsoft\Windows\CurrentVersion\App Paths\uninst.exe"
!define PRODUCT_UNINST_KEY "Software\Microsoft\Windows\CurrentVersion\Uninstall\${PRODUCT_NAME}"
!define PRODUCT_UNINST_ROOT_KEY "HKLM"
!define PRODUCT_STARTMENU_REGVAL "NSIS:StartMenuDir"

SetCompressor lzma

; MUI 1.67 compatible ------
!include "MUI.nsh"

; MUI Settings
!define MUI_ABORTWARNING
!define MUI_ICON "${NSISDIR}\Contrib\Graphics\Icons\modern-install-colorful.ico"
!define MUI_UNICON "${NSISDIR}\Contrib\Graphics\Icons\modern-uninstall-colorful.ico"

; Some variables
Var ICONS_GROUP
Var STARTDIR

; Welcome page (allow 3 lines if long pre-release)
!define MUI_WELCOMEPAGE_TITLE_3LINES
!define MUI_WELCOMEFINISHPAGE_BITMAP_NOSTRETCH
!define MUI_WELCOMEFINISHPAGE_BITMAP "C:\msys64\home\bernhard\installer\coot_pic_welcome.bmp"
!define MUI_HEADERIMAGE
!define MUI_HEADERIMAGE_BITMAP_NOSTRETCH
!define MUI_HEADERIMAGE_BITMAP "C:\msys64\home\bernhard\installer\coot_pic_header.bmp"
!define MUI_PAGE_CUSTOMFUNCTION_PRE InitPreFunction
!insertmacro MUI_PAGE_WELCOME
; Components page
!define MUI_PAGE_CUSTOMFUNCTION_PRE CheckForUpdate
!insertmacro MUI_PAGE_COMPONENTS
; License page
!define MUI_PAGE_CUSTOMFUNCTION_PRE CheckForUpdate
!define MUI_LICENSEPAGE_TEXT_TOP "Coot licence"
!define MUI_LICENSEPAGE_CHECKBOX
!insertmacro MUI_PAGE_LICENSE "C:\msys64\home\bernhard\autobuild\extras\COPYING"
!define MUI_PAGE_CUSTOMFUNCTION_PRE SkipProbeReduceInfo
!define MUI_LICENSEPAGE_TEXT_TOP "Probe && Reduce information"
!define MUI_LICENSEPAGE_TEXT_BOTTOM "This is just for your information. No licence as such is required."
!define MUI_LICENSEPAGE_BUTTON "Got it!"
!insertmacro MUI_PAGE_LICENSE "C:\msys64\home\bernhard\autobuild\extras\probe_reduce.txt"
; Directory page (install directory)
!define MUI_PAGE_CUSTOMFUNCTION_PRE ComponentPost
!define MUI_PAGE_CUSTOMFUNCTION_SHOW DirectoryShow
!insertmacro MUI_PAGE_DIRECTORY
; Start menu page
!define MUI_STARTMENUPAGE_NODISABLE
!define MUI_STARTMENUPAGE_DEFAULTFOLDER "WinCoot"
!define MUI_STARTMENUPAGE_REGISTRY_ROOT "${PRODUCT_UNINST_ROOT_KEY}"
!define MUI_STARTMENUPAGE_REGISTRY_KEY "${PRODUCT_UNINST_KEY}"
!define MUI_STARTMENUPAGE_REGISTRY_VALUENAME "${PRODUCT_STARTMENU_REGVAL}"
!define MUI_PAGE_CUSTOMFUNCTION_PRE CheckForUpdate
!insertmacro MUI_PAGE_STARTMENU Application $ICONS_GROUP
; Directory page to set the start directory
!define MUI_DIRECTORYPAGE_VARIABLE $STARTDIR
!define MUI_DIRECTORYPAGE_TEXT_DESTINATION $STARTDIR
!define MUI_PAGE_CUSTOMFUNCTION_SHOW DirectoryShow
!define MUI_PAGE_CUSTOMFUNCTION_PRE CheckForUpdate
!insertmacro MUI_PAGE_DIRECTORY

; Instfiles page
!insertmacro MUI_PAGE_INSTFILES
; Finish page
; run the bat file changes and gdk-pixbuf before we finish the installation
!define MUI_PAGE_CUSTOMFUNCTION_PRE FinishPagePreFunction
!define MUI_FINISHPAGE_RUN "$INSTDIR\wincoot.bat"
;!define MUI_FINISHPAGE_RUN_PARAMETERS '/c "$INSTDIR\wincoot.bat"'
!define MUI_FINISHPAGE_TITLE_3LINES
!insertmacro MUI_PAGE_FINISH

; Uninstaller pages
!insertmacro MUI_UNPAGE_INSTFILES

; Language files
!insertmacro MUI_LANGUAGE "English_WinCoot"
;!include "${NSISDIR}\Contrib\Language Files\English_WinCoot.nsh"

; for replacing strings
!include "WordFunc.nsh"
!insertmacro WordReplace
; actually a macro wrapping the error handler function
; well, doesnt seem to work and always throws errors
; itself somewhere, so revert to old school if and repeating
; code. Keep in case it can be resolved at some point...
!macro MyErrorHandlerMacro un

Function ${un}ErrorHandlerFunction

  Pop $0  ; Abort Flag
  Pop $1  ; Message txt
  Pop $2  ; error code

; convert strings to ints
  IntOp $0 $0 + 0
  IntOp $2 $2 + 0

; the message box could be rather a question to abort!?
  IfSilent +2 0
  MessageBox MB_OK 'Error in Installation. Aborting!$\n$\r$\n$\r$1'

  DetailPrint $1
  SetErrorLevel $2
  ${If} $0 == 1
    Abort $1
  ${EndIf}

FunctionEnd
!macroend

; for error handling, need the second one for uninstall section
!define ErrorHandler "!insertmacro ErrorHandler"
!define UnErrorHandler "!insertmacro UnErrorHandler"

!insertmacro MyErrorHandlerMacro ""
!insertmacro MyErrorHandlerMacro "un."


; ErrorCode is returned by the installer, ErrorText is written in log
; AbortFlag is 1/"True" or 0/"False", (un)installer terminates upon error
!macro ErrorHandler ErrorCode ErrorText AbortFlag
  Push "${ErrorCode}"
  Push "${ErrorText}"
  Push "${AbortFlag}"
  Call ErrorHandlerFunction
!macroend

; same for uninstall (needs "un." in front of function)
!macro UnErrorHandler ErrorCode ErrorText AbortFlag
  Push "${ErrorCode}"
  Push "${ErrorText}"
  Push "${AbortFlag}"
  Call un.ErrorHandlerFunction
!macroend

; MUI end ------

Name "${PRODUCT_NAME} ${PRODUCT_VERSION}"

!ifndef binary_dir
; default to pre-release output dir
!define binary_dir "C:\msys64\home\bernhard\public_html\software\binaries\pre-release"
!endif
; just in case something goes wrong, we define src to be pre-release (default)
!ifndef src_dir
!define src_dir "C:\msys64\home\bernhard\autobuild\MINGW64_NT-6.1-bernie-pre-release-gtk2-shared"
!endif
!define outputname "${binary_dir}\${PRODUCT_NAME}-${PRODUCT_VERSION}.exe"
OutFile "${outputname}"

XPStyle on

InstallDir "C:\WinCoot"
InstallDirRegKey HKLM "${PRODUCT_DIR_REGKEY}" ""
ShowInstDetails show
ShowUnInstDetails show



#####################
# SECTIONS
#####################

Section "!WinCoot" SEC01
  ClearErrors
  SectionIn RO
  ; first check for Admin rights and if we shall install for all users
  ; change logic here. Install for everyone if possible otherwise only
  ; for current user. Done later in Shortcut section. Based on LilyPond
  ;Var /GLOBAL install_all_users
  ;StrCpy $install_all_users "False"
  ;userInfo::getAccountType
  ;pop $0
  ;StrCmp $0 "Admin" +2
  ;Goto install_single_user
  ; if silent installer, install for all
  ;IfSilent install_all_users
  ;messageBox MB_ICONQUESTION|MB_YESNO "Detected Administrator rights.$\nDo you want to install WinCoot for all Users?" IDYES install_all_users IDNO install_single_user
  ;install_all_users:
  ;  SetShellVarContext all
  ;  StrCpy $install_all_users "True"
  ;  Goto continue_user
  ;install_single_user:
  ;  SetShellVarContext current
  ;continue_user:

;  !insertmacro BIMAGE "C:\msys64\home\bernhard\installer\coot_pic.bmp" /RESIZETOFIT
  SetOverwrite ifnewer
  SetOutPath "$INSTDIR"

; check if wincoot.bat exists, if so ask if overwrite (or not)
; New logic:
; check if wincoot.bat exists, if so make a backup copy
; (this is a tmp file, doesnt seem to be deleted, maybe virus thing? not sure what's happening)
; install new wincoot.bat
; deal with rest in FinishPagePreFunction
  Var /GLOBAL have_bat
  StrCpy $have_bat "False"
  IfFileExists "$INSTDIR\runwincoot.bat" "" endifoldbat
    Rename $INSTDIR\runwincoot.bat $INSTDIR\wincoot.bat
  endifoldbat:

  IfFileExists "$INSTDIR\wincoot.bat" "" endifbat
    StrCpy $have_bat "True"
  endifbat:

  SetOverwrite on
  File /oname=$INSTDIR\wincoot.bat.tmp "C:\msys64\home\bernhard\Projects\coot\windows\wincoot.bat"
  File /oname=$INSTDIR\wincoot-for-ccp4i2.bat.tmp "C:\msys64\home\bernhard\Projects\coot\windows\wincoot-for-ccp4i2.bat"

  SetOverwrite ifnewer
; bin DIR
  SetOutPath "$INSTDIR\bin"
  SetOverwrite on
  File "${src_dir}\bin\coot-bin.exe"
  File "${src_dir}\bin\coot-density-score-by-residue-bin.exe"
  File "${src_dir}\bin\findligand-bin.exe"
  File "${src_dir}\bin\findwaters-bin.exe"
  File "${src_dir}\bin\mini-rsr-bin.exe"
  File "${src_dir}\bin\dynarama-bin.exe"
  File "${src_dir}\bin\pyrogen.bat"
  File "${src_dir}\bin\python.exe"
  File "C:\msys64\home\bernhard\Projects\coot\windows\dynarama.bat"
  File "C:\msys64\home\bernhard\Projects\coot\windows\findligand.bat"
  File "C:\msys64\home\bernhard\Projects\coot\windows\findwaters.bat"
  SetOverwrite ifnewer
  File "C:\msys64\home\bernhard\autobuild\extras\coot-icon.ico"
  File "C:\msys64\home\bernhard\autobuild\extras\rama_all.ico"
  File "${src_dir}\bin\*.dll"
  File "${src_dir}\bin\coot-bfactan.exe"
  File "${src_dir}\bin\coot"
  File "${src_dir}\bin\coot-mini-rsr"
  File "${src_dir}\bin\coot-available-comp-id.exe"
  File "${src_dir}\bin\coot-compare-dictionaries.exe"
;  File "${src_dir}\bin\coot-dictionary-bond-distributions.exe"
  File "${src_dir}\bin\coot-make-shelx-restraints.exe"
  File "${src_dir}\bin\coot-density-score-by-residue"
  File "${src_dir}\bin\findligand"
  File "${src_dir}\bin\findwaters"
  File "${src_dir}\bin\coot-fix-nomenclature-errors.exe"
  File "${src_dir}\bin\gdk-pixbuf-query-loaders.exe"
  ; render (or more?) from raster3d?!
  File "${src_dir}\bin\render.exe"
  ;gunzip needed?? dont think so
  File "C:\msys64\home\bernhard\autobuild\extras\gunzip"
  File "C:\msys64\home\bernhard\autobuild\extras\gzip.exe"
  File "${src_dir}\bin\lidia.exe"
  ;still needed?
;  File "C:\msys64\home\bernhard\autobuild\extras\msvcr90.dll"
;  File "C:\msys64\home\bernhard\autobuild\extras\Microsoft.VC90.CRT.manifest"
;  File "${src_dir}\bin\pkg-config.exe"
  File "C:\msys64\home\bernhard\autobuild\extras\ppm2bmp.exe"
  ; now the new mingw files for static compilation (or not)
  ; FIXME:: seems to be only needed in newer versions of msys
  ; ....... and maybe if we have shared compilation - not yet
  ; 19/2/19 take the system file on the system its build on
  ; included in bundling
  ;File "C:\MinGW\bin\libstdc++-6.dll"
  ;File "C:\MinGW\bin\libgcc_s_dw2-1.dll"
; PYTHON stuff new
; now in lib/python2.7 see below
  ;SetOutPath "$INSTDIR\python27"
  ;File /r "${src_dir}\python27\*.*"
; etc things - not all needed any more?!?!? FIXME
;  SetOutPath "$INSTDIR\etc\fonts"
;  File "${src_dir}\etc\fonts\*"
;  SetOutPath "$INSTDIR\etc\gtk-2.0"
;  File "${src_dir}\etc\gtk-2.0\gtk.immodules"
;  SetOutPath "$INSTDIR\etc\pango"
;  File "${src_dir}\etc\pango\pango.modules"
; for ccp4 version
  SetOutPath "$INSTDIR\etc"
  File "C:\msys64\home\bernhard\Projects\coot\windows\runwincoot_ccp4.bat"
  File "C:\msys64\home\bernhard\Projects\coot\windows\runwincoot_ccp4_vista.bat"
; SHARE
  SetOutPath "$INSTDIR\share\icons"
  File /r "${src_dir}\share\icons\*.*"
  SetOutPath "$INSTDIR\share\coot"
  File /r "${src_dir}\share\coot\*.*"
  ;lib
  ; maybe the boost and rdkit dlls should be in bin rather than lib?!
  SetOutPath "$INSTDIR\lib"
  File /r "${src_dir}\lib\*.dll"
  SetOutPath "$INSTDIR\lib\python2.7"
  File /r "${src_dir}\lib\python2.7\*.*"
  SetOutPath "$INSTDIR\lib\gdk-pixbuf-2.0"
  File /r "${src_dir}\lib\gdk-pixbuf-2.0\*.*"
  SetOutPath "$INSTDIR\lib\gtk-2.0"
  File /r "${src_dir}\lib\gtk-2.0\*.*"
!ifdef WITH_GUILE
  ;guile things (shouldnt they be in Guile section?! Never mind)
  SetOutPath "$INSTDIR\bin"
  File "${src_dir}-guile\bin\*guile*"
  File "${src_dir}-guile\bin\libgmp-3.dll"
  File "${src_dir}-guile\bin\readline5.dll"
  File "${src_dir}-guile\bin\libltdl3.dll"
  SetOutPath "$INSTDIR\share\guile\1.8"
  File "${src_dir}-guile\share\guile\1.8\*"
  SetOutPath "$INSTDIR\share\guile\1.8\ice-9"
  File "${src_dir}-guile\share\guile\1.8\ice-9\*"
  SetOutPath "$INSTDIR\share\guile\1.8\ice-9\debugger"
  File "${src_dir}-guile\share\guile\1.8\ice-9\debugger\*"
  SetOutPath "$INSTDIR\share\guile\1.8\lang\elisp"
  File "${src_dir}-guile\share\guile\1.8\lang\elisp\*"
  SetOutPath "$INSTDIR\share\guile\1.8\lang\elisp\internals"
  File "${src_dir}-guile\share\guile\1.8\lang\elisp\internals\*"
  SetOutPath "$INSTDIR\share\guile\1.8\lang\elisp\primitives"
  File "${src_dir}-guile\share\guile\1.8\lang\elisp\primitives\*"
  SetOutPath "$INSTDIR\share\guile\1.8\oop"
  File "${src_dir}-guile\share\guile\1.8\oop\*"
  SetOutPath "$INSTDIR\share\guile\1.8\oop\goops"
  File "${src_dir}-guile\share\guile\1.8\oop\goops\*"
  SetOutPath "$INSTDIR\share\guile\1.8\scripts"
  File "${src_dir}-guile\share\guile\1.8\scripts\*"
  SetOutPath "$INSTDIR\share\guile\1.8\srfi"
  File "${src_dir}-guile\share\guile\1.8\srfi\*"
  SetOutPath "$INSTDIR\share\guile\gtk"
  File "${src_dir}-guile\share\guile\gtk\*"
  SetOutPath "$INSTDIR\share\guile\gtk-2.0"
  File "${src_dir}-guile\share\guile\gtk-2.0\*"
  SetOutPath "$INSTDIR\share\guile\gui"
  File "${src_dir}-guile\share\guile\gui\*"
  SetOutPath "$INSTDIR\share\guile\site"
  File "${src_dir}-guile\share\guile\site\*"
!endif
; WITH_GUILE
  ;docs
  SetOutPath "$INSTDIR\doc"
  File "C:\msys64\home\bernhard\autobuild\extras\coot-user-manual.pdf"
  File "C:\msys64\home\bernhard\autobuild\extras\crib-sheet.pdf"
  File "C:\msys64\home\bernhard\autobuild\extras\tutorial.pdf"
  File "C:\msys64\home\bernhard\autobuild\extras\tutorial-2.pdf"
  ;;secondary structure(s)
  ;SetOutPath "$INSTDIR\share\coot\ss-reference-structures"
  ;File "C:\msys64\home\bernhard\autobuild\extras\ss-reference-structures\*"
  ;make a backupdir, so that COOT_BACKUP_DIR has a defined directory
  SetOutPath "$INSTDIR\coot-backup"
  ; set outpath to $INSTDIR so that shortcuts are started in $INSTDIR
  SetOutPath "$INSTDIR"
  IfErrors 0 +6
;    ${ErrorHandler} 1 "Error in installation. Could not write files." 1
     DetailPrint "Error in installation. Could not write files."
     SetErrorLevel 1
     IfSilent +2 0
       MessageBox MB_OK 'Error in Installation. Aborting!$\n$\r$\n$\rCould not write files.'
     Abort "Abort. Could not write files."

SectionEnd


Section /o "Windows feel" SEC02
  ClearErrors
  SetOverwrite on
  SetOutPath "$INSTDIR\share\coot"
  File "C:\msys64\home\bernhard\autobuild\extras\cootrc"
  SetOverwrite ifnewer
;  maybe here the other guile things?!
  IfErrors 0 +5
  ;  ${ErrorHandler} 2 "Error in installation. Could not install Windows feel." 1
     DetailPrint "Error in installation. Could not install Windows feel. Continuing."
     SetErrorLevel 2
     IfSilent +2 0
        MessageBox MB_OK 'Error in Installation. Continuing!$\n$\r$\n$\rCould not install Windows feel'

SectionEnd

; probe and reduce here
Section /o "Add probe&reduce" SEC03
  ClearErrors
  SetOverwrite on
  ; This needs to go in another bin as to void conflicst with the libstd and
  ; libgcc dlls. Needs addition to PATH as well
  SetOutPath "$INSTDIR\bin\extras"
  File "C:\msys64\home\bernhard\autobuild\extras\probe.exe"
  File "C:\msys64\home\bernhard\autobuild\extras\reduce.exe"
  File "C:\msys64\home\bernhard\autobuild\extras\libstdc++-6.dll"
  File "C:\msys64\home\bernhard\autobuild\extras\libgcc_s_dw2-1.dll"
  SetOutPath "$INSTDIR\share\coot"
  File "C:\msys64\home\bernhard\autobuild\extras\reduce_wwPDB_het_dict.txt"
  SetOverwrite ifnewer

  IfErrors 0 +5
    ; ${ErrorHandler} 2 "Error in installation. Could not install probe&reduce." 1
     DetailPrint "Error in installation. Could not install probe&reduce. Continuing."
     SetErrorLevel 2
     IfSilent +2 0
        MessageBox MB_OK 'Error in Installation. Continuing!$\n$\r$\n$\rCould not install probe & reduce'

SectionEnd

!ifdef WITH_GUILE
; we dont want guile for now
; 3rd section for guile
Section /o "Guile/Scheme Add-On" SEC04
  SetOverwrite on
  SetOutPath "$INSTDIR\bin"
  File "${src_dir}-guile\bin\coot-real.exe"
  SetOverwrite ifnewer
  maybe here the other guile things?!
SectionEnd
!endif
; WITH_GUILE

Section -AddIcons
  ;; First install for all users, if anything fails, install
  ;; for current user only.
  ClearErrors

  SetShellVarContext all
  ; let's see what happens when we try for all
  Call MakeIcons
  ClearErrors

  ; if error delete what may be there (but there shouldnt be anything?)
  IfErrors 0 exit
  SetShellVarContext current
  Call MakeIcons

 exit:
   SetShellVarContext current
  ; in case we want to install silently (no user intervention)
  ; we have to call the finish page function, otherwise runwincoot.bat
  ; won't be edited.
  IfSilent 0 +2
    Call FinishPagePreFunction
  IfErrors 0 +5
    ; ${ErrorHandler} 3 "Error in installation. Could not install icons." 0
    DetailPrint "Error in installation. Could not install icons. Continuing."
    SetErrorLevel 3
    IfSilent +2 0
        MessageBox MB_OK 'Error in Installation. Continuing!$\n$\r$\n$\rCould not install icons.'
SectionEnd

Section -Post
  ClearErrors
  WriteUninstaller "$INSTDIR\uninst.exe"
;  NO MESSING WITH THE REGISTRY!!!!
;  WriteRegStr HKLM "${PRODUCT_DIR_REGKEY}" "" "$INSTDIR\uninst.exe"
;  WriteRegStr ${PRODUCT_UNINST_ROOT_KEY} "${PRODUCT_UNINST_KEY}" "DisplayName" "$(^Name)"
;  WriteRegStr ${PRODUCT_UNINST_ROOT_KEY} "${PRODUCT_UNINST_KEY}" "UninstallString" "$INSTDIR\uninst.exe"
;  WriteRegStr ${PRODUCT_UNINST_ROOT_KEY} "${PRODUCT_UNINST_KEY}" "DisplayIcon" "$INSTDIR\uninst.exe"
;  WriteRegStr ${PRODUCT_UNINST_ROOT_KEY} "${PRODUCT_UNINST_KEY}" "DisplayVersion" "${PRODUCT_VERSION}"
;  WriteRegStr ${PRODUCT_UNINST_ROOT_KEY} "${PRODUCT_UNINST_KEY}" "URLInfoAbout" "${PRODUCT_WEB_SITE}"
;  WriteRegStr ${PRODUCT_UNINST_ROOT_KEY} "${PRODUCT_UNINST_KEY}" "Publisher" "${PRODUCT_PUBLISHER}"
  IfErrors 0 +5
    ; ${ErrorHandler} 4 "Error in installation. Could not write uninstaller." 0
    DetailPrint "Error in installation. Could not write uninstaller."
    SetErrorLevel 4
    IfSilent +2 0
        MessageBox MB_OK 'Error in Installation. Continuing!$\n$\r$\n$\rCould not write uninstaller.'
SectionEnd



; Section descriptions
!insertmacro MUI_FUNCTION_DESCRIPTION_BEGIN
  !insertmacro MUI_DESCRIPTION_TEXT ${SEC01} "This is 'default' WinCoot (${WinCootVersion}) $\n$\nPython scripting only"
  !insertmacro MUI_DESCRIPTION_TEXT ${SEC02} "Tick if you want a $\nWindowsy feeling to WinCoot"
  !insertmacro MUI_DESCRIPTION_TEXT ${SEC03} "Select if you want $\nprobe and reduce installed.$\nNote: Usually not required if you have CCP4 installed."
; disable guile for now
!ifdef WITH_GUILE
  !insertmacro MUI_DESCRIPTION_TEXT ${SEC04} "Tick if you want additionally $\nGuile/Scheme scripting"
!endif
!insertmacro MUI_FUNCTION_DESCRIPTION_END


Section Uninstall
  ClearErrors
  Delete "$INSTDIR\${PRODUCT_NAME}.url"
  Delete "$INSTDIR\uninst.exe"
  RMDir /r "$INSTDIR\share"
  ; keep the next 2 in case it was there from previous installations
  Delete "$INSTDIR\bin\Lib\*"
  Delete "$INSTDIR\bin\DLLs\*"
  ; only remove installed files in bin
  ; (in case there exists probe, reduce e.g.)
  Delete "$INSTDIR\bin\coot-icon.ico"
  Delete "$INSTDIR\bin\*.dll"
  Delete "$INSTDIR\bin\bfactan.exe"
  Delete "$INSTDIR\bin\coot-bfactan.exe"
  Delete "$INSTDIR\bin\coot"
  Delete "$INSTDIR\bin\coot-real.exe"
  Delete "$INSTDIR\bin\coot-bin.exe"
  Delete "$INSTDIR\libexec\*.exe"
  Delete "$INSTDIR\bin\coot-density-score-by-residue"
  Delete "$INSTDIR\bin\density-score-by-residue-bin.exe"
  Delete "$INSTDIR\bin\density-score-by-residue"
  Delete "$INSTDIR\bin\dynarama-bin.exe"
  Delete "$INSTDIR\bin\dynarama.bat"
  Delete "$INSTDIR\bin\dynarama"
  Delete "$INSTDIR\bin\pyrogen.bat"
  Delete "$INSTDIR\bin\python.exe"
  Delete "$INSTDIR\bin\rama_all.ico"
  Delete "$INSTDIR\bin\findligand"
  Delete "$INSTDIR\bin\findligand-bin.exe"
  Delete "$INSTDIR\bin\findligand-real.exe"
  Delete "$INSTDIR\bin\findligand.bat"
  Delete "$INSTDIR\bin\findwaters"
  Delete "$INSTDIR\bin\findwaters-bin.exe"
  Delete "$INSTDIR\bin\findwaters-real.exe"
  Delete "$INSTDIR\bin\findwaters.bat"
  Delete "$INSTDIR\bin\coot-available-comp-id.exe"
  Delete "$INSTDIR\bin\coot-compare-dictionaries.exe"
  Delete "$INSTDIR\bin\coot-density-score-by-residue-bin.exe"
  Delete "$INSTDIR\bin\coot-fix-nomenclature-errors.exe"
  Delete "$INSTDIR\bin\coot-make-shelx-restraints.exe"
  Delete "$INSTDIR\bin\coot-mini-rsr"
  Delete "$INSTDIR\bin\fix-nomenclature-errors.exe"
  Delete "$INSTDIR\bin\gdk-pixbuf-csource.exe"
  Delete "$INSTDIR\bin\gdk-pixbuf-query-loaders.exe"
  Delete "$INSTDIR\bin\glib-genmarshal.exe"
  Delete "$INSTDIR\bin\glib-gettextize"
  Delete "$INSTDIR\bin\glib-mkenums"
  Delete "$INSTDIR\bin\gobject-query.exe"
  Delete "$INSTDIR\bin\gsl-config"
  Delete "$INSTDIR\bin\gspawn-win32-helper-console.exe"
  Delete "$INSTDIR\bin\gspawn-win32-helper.exe"
  Delete "$INSTDIR\bin\gtk-builder-convert"
  Delete "$INSTDIR\bin\gtk-demo.exe"
  Delete "$INSTDIR\bin\gtk-query-immodules-2.0.exe"
  Delete "$INSTDIR\bin\gunzip"
  Delete "$INSTDIR\bin\gzip.exe"
  Delete "$INSTDIR\bin\iconv.exe"
  Delete "$INSTDIR\bin\lidia.exe"
  Delete "$INSTDIR\bin\mini-rsr-bin.exe"
  Delete "$INSTDIR\bin\render.exe"
  Delete "$INSTDIR\bin\Microsoft.VC90.CRT.manifest"
  Delete "$INSTDIR\bin\pango-querymodules.exe"
  Delete "$INSTDIR\bin\pango-view.exe"
  Delete "$INSTDIR\bin\pkg-config.exe"
  Delete "$INSTDIR\bin\ppm2bmp.exe"
  ; probe & reduce (optional?)
  Delete "$INSTDIR\bin\extras\*"
  ; guile things
  Delete "$INSTDIR\bin\*guile*"
  Delete "$INSTDIR\bin\libgmp-3.dll"
  Delete "$INSTDIR\bin\readline5.dll"
  Delete "$INSTDIR\bin\libltdl3.dll"
  ; continue deleting everything
  RMDir /r "$INSTDIR\etc"
  RMDir /r "$INSTDIR\lib"
  RMDir /r "$INSTDIR\examples"
  RMDir /r "$INSTDIR\doc"
  ; we shall only remove installed files here (to e.g. keep .coot.py untouched)
  Delete "$INSTDIR\*wincoot*.bat"
  Delete "$INSTDIR\Coot.url"

  !insertmacro MUI_STARTMENU_GETFOLDER "Application" $ICONS_GROUP
  SetShellVarContext all
  Delete "$SMPROGRAMS\$ICONS_GROUP\*"
  Delete "$DESKTOP\WinCoot.lnk"
  Delete "$QUICKLAUNCH\WinCoot.lnk"
  RMDir "$SMPROGRAMS\$ICONS_GROUP"
  SetShellVarContext current
  Delete "$SMPROGRAMS\$ICONS_GROUP\*"
  Delete "$DESKTOP\WinCoot.lnk"
  Delete "$QUICKLAUNCH\WinCoot.lnk"
  RMDir "$SMPROGRAMS\$ICONS_GROUP"

  ; keep these in case they exist from previous versions
  RMDir /r "$INSTDIR\bin\Lib"
  RMDir "$INSTDIR\bin\DLLs"
  RMDir "$INSTDIR\bin\extras"
  RMDir "$INSTDIR\bin"
  RMDir /r "$INSTDIR\libexec"
  ; FIXME:: what to do with python?
  ; this is only for old installs...
  RMDir /r "$INSTDIR\python27"
  ; will only be removed if empty
  RMDir "$INSTDIR\.coot-preferences"
  RMDir "$INSTDIR\coot-backup"
  RMDir "$INSTDIR\"

;  DeleteRegKey ${PRODUCT_UNINST_ROOT_KEY} "${PRODUCT_UNINST_KEY}"
;  DeleteRegKey HKLM "${PRODUCT_DIR_REGKEY}"
  SetAutoClose true
  IfErrors 0 +5
    ;${UnErrorHandler} 5 "Error in uninstallation. Could not completely uninstall." 1
    DetailPrint "Error in uninstallation. Could not completely uninstall."
    SetErrorLevel 5
    IfSilent +2 0
        MessageBox MB_OK 'Error in Uninstallation. Aborting!$\n$\r$\n$\rCould not completely uninstall.'
SectionEnd

######################
# 'BUILD-IN' FUNCTIONS
######################

# BL says:: disable for now, since no guile available anyway
!ifdef WITH_GUILE
Function .onSelChange
   ${If} ${SectionIsSelected} ${SEC03}
      MessageBox MB_OK|MB_ICONEXCLAMATION "You have with guile selected. Sure? This may not work perfectly.$\r$\n\"
   ${EndIf}
FunctionEnd
!endif

Function .onInit
  ClearErrors

  ; logging
  ; if old log exists save simply (could be by date)
  ; StrCpy $TEMP .
  ; Dont use INSTDIR as it doesnt exist from the beginning use $TEMP instead.
  ; Maybe use a different place at some point.
  IfFileExists $TEMP\install.log 0 +4
    IfFileExists $TEMP\install.log.1 0 +2
      Delete $TEMP\install.log.1
    Rename $TEMP\install.log $TEMP\install.log.1
  ; Seems install.log always goes to $INSTDIR, since this doesnt
  ; exists yet we try to set it in another way
  Push $INSTDIR
  StrCpy $INSTDIR $TEMP
  LogSet on
  Pop $INSTDIR

    ; Get Command line parameters
        var /GLOBAL INSTDIR_TMP
        StrCpy $INSTDIR_TMP $INSTDIR
        StrCpy $INSTDIR ""

	var /GLOBAL cmdLineParams
	Push $R0

	${GetParameters} $cmdLineParams

	; /? param (help)
	ClearErrors
	${GetOptions} $cmdLineParams '/?' $R0
	IfErrors +3 0
        Call ShowOptions
	Abort

	ClearErrors
	${GetOptions} $cmdLineParams '/help' $R0
	IfErrors +3 0
        Call ShowOptions
	Abort

	Pop $R0

	ClearErrors

    ; Initialise options

	Var /GLOBAL update
	Var /GLOBAL delete_installer
	Var /GLOBAL start_coot

	StrCpy $update             	  0
	StrCpy $delete_installer      	  0
	StrCpy $start_coot      	  0

    ; Parse Parameters

	Push $R0
	Call parseParameters
	Pop $R0

        ${If} $INSTDIR == ""
          ${If} $update = 1
           Call NoInstDirGiven
          ${Else}
           StrCpy $INSTDIR $INSTDIR_TMP
          ${EndIf}
        ${EndIf}
; insert params END
; no more examples dir (have data dir now)
; so default start dir is $INSTDIR
  ${If} $STARTDIR == ""
    StrCpy $STARTDIR "$INSTDIR"
  ${EndIf}
  IfErrors 0 +6
    ;${ErrorHandler} 6 "Error in installation. Could not initiate installation." 1
    DetailPrint "Error in installation. Could not initiate installation."
    SetErrorLevel 6
    IfSilent +2 0
        MessageBox MB_OK 'Error in Installation. Aborting!$\n$\r$\n$\rCould not initiate installation.'
    Abort "Error in installation. Could not initiate installation. Aborting."

FunctionEnd

; finallising
Function .onGUIEnd
  ClearErrors
  ; remove tmp bat file if exists (Delete wont do, need to do it via DOS shell command)
  ; possibly an anti virus thing?!?
  ; delete both bat files i.e. for i2 too
  IfFileExists "$INSTDIR\wincoot.bat.tmp" "" cont
     ;Delete "$INSTDIR\wincoot.bat.tmp"
     Exec 'cmd /c del "$INSTDIR\wincoot.bat.tmp"'
  cont:
  IfFileExists "$INSTDIR\wincoot-for-ccp4i2.bat.tmp" "" conti2
     Exec 'cmd /c del "$INSTDIR\wincoot-for-ccp4i2.bat.tmp"'
  conti2:

  ; delete the installer
  ${If} $delete_installer = 1
   !insertmacro StrStr $0 "$EXEDIR" "pending-install"
   ${If} $0 == ""
    Exec 'cmd /c del "$EXEDIR\${PRODUCT_NAME}-${PRODUCT_VERSION}.exe"'
   ${Else}
    Exec 'cmd /c rmdir /s /q "$EXEDIR"'
   ${EndIf}
  ${EndIf}

  ; Finally if run then start WinCoot in StartDir
  ${If} $start_coot = 1
    SetOutPath $STARTDIR
    Exec $INSTDIR\wincoot.bat
  ${EndIf}
  IfErrors 0 +6
    ; ${ErrorHandler} 7 "Error in installation. Could not write/edit runwincoot.bat." 1
    DetailPrint "Error in installation. Could not write/edit runwincoot.bat."
    SetErrorLevel 7
    IfSilent +2 0
        MessageBox MB_OK 'Error in Installation. Aborting!$\n$\r$\n$\rCould not write/edit runwincoot.bat.'
    Abort "Error in installation. Could not write/edit (run)wincoot.bat. Aborting!"

FunctionEnd

Function un.onUninstSuccess
  HideWindow
  MessageBox MB_ICONINFORMATION|MB_OK "$(^Name) was successfully removed from your computer." /SD IDOK
FunctionEnd

Function un.onInit
  MessageBox MB_ICONQUESTION|MB_YESNO|MB_DEFBUTTON2 "Are you sure you want to completely remove $(^Name) and all of its components?" /SD IDYES IDYES +2
  Abort
FunctionEnd

###################
# GENERAL FUNCTIONS
###################

Function MakeBackground
  ; background depends on update
  ${If} $update = 0
    ; from dark green to dark red : 0 0x80 0 0x80 0 0
    BgImage::SetBg /GRADIENT 0 0x80 0 0x80 0 0
    BgImage::AddText "$(^Name) Installer"
    ;BgImage::AddImage C:\msys64\home\bernhard\installer\coot_pic_welcome.bmp 150 0
    CreateFont $R0 "Comic Sans MS" 30 700
    # add a blue shadow for the text
    BgImage::AddText "$(^Name) Installer" $R0 0 0 255 48 48 798 198
    # add a green shadow for the text
    BgImage::AddText "$(^Name) Installer" $R0 0 255 0 52 52 802 202
    # add the text
    BgImage::AddText "$(^Name) Installer" $R0 255 0 0 50 50 800 200
    BgImage::Redraw
  ${EndIf}
FunctionEnd

Function parseParameters

    ; /instdir=
    ; What about spaces in dir name?? USE quotes!?
    ${GetOptions} $cmdLineParams '/instdir=' $R0
    IfErrors continue0 0
    ${If} $R0 == ""
    MessageBox MB_ICONSTOP "No install directory given.$\n$\n\
               Please use command line option /instdir=instdir$\n$\n\
               Abort installation!"
    Abort
    ${Else}
    StrCpy $INSTDIR $R0
    ${EndIf}
    continue0:

    ; /update
    ${GetOptions} $cmdLineParams '/update' $R0
    IfErrors continue1 0
    StrCpy $update 1
    continue1:

    ; /autoupdate
    ${GetOptions} $cmdLineParams '/autoupdate' $R0
    IfErrors continue2 0
    StrCpy $update 1
    StrCpy $delete_installer 1
    StrCpy $start_coot 1
    continue2:

    ; /startdir=
    ; What about spaces in dir name?? USE quotes!?
    ${GetOptions} $cmdLineParams '/startdir=' $R0
    IfErrors continue3 0
    ${If} $R0 == ""
    MessageBox MB_ICONSTOP "No start directory given.$\n$\n\
               Please use command line option /startdir=startdir$\n$\n\
               Abort installation!"
    Abort
    ${Else}
    StrCpy $STARTDIR $R0
    ${EndIf}
    continue3:

    ; /run
    ${GetOptions} $cmdLineParams '/run' $R0
    IfErrors +2 0
    StrCpy $start_coot 1

    ; /del_inst
    ${GetOptions} $cmdLineParams '/del_inst' $R0
    IfErrors +2 0
    StrCpy $delete_installer 1

FunctionEnd

Function NoInstDirGiven
    IfSilent 0 +2
    MessageBox MB_ICONSTOP "No installation directory given.$\n$\n\
               Please use command line option /instdir=instdir$\n$\n\
               Abort installation!"
    Abort
FunctionEnd

Function ShowOptions
	MessageBox MB_ICONQUESTION "WinCoot Installer Command Line Help$\n$\n\
        Options:$\n\
          /D=instdir         - defines installation directory to instdir$\n\
          /instdir=instdir   - as above (preferred!!! needs to be first)$\n\
          /startdir=statdir  - defines dir where Coot will start (after install/update)$\n\
          /update            - update WinCoot (no questions asked, requires /D)$\n\
          /autoupdate        - as /update but removes installer and starts WinCoot$\n\
          /del_inst          - removes the installer after installation$\n\
          /run               - run WinCoot after installation$\n\
          /S                 - silent installation (and uninstall)$\n\
          /help | /?         - this information"
FunctionEnd

; Dont show pages when updating!
Function CheckForUpdate
  ${If} $update = 1
     Abort
  ${EndIf}
FunctionEnd

; to skip the probe & reduce info if we dont install this
Function SkipProbeReduceInfo
  ${IfNot} ${SectionIsSelected} ${SEC03}
     abort ; no selected
  ${EndIf}
FunctionEnd

; before we start pages
Function InitPreFunction
  Call CheckForUpdate
  Call MakeBackground
FunctionEnd

; inistiallize the directory page count
Function ComponentPost
  Call CheckForUpdate
  Var /GLOBAL DIRCOUNT
  StrCpy $DIRCOUNT "0"
FunctionEnd

Function DirectoryShow
  StrCmp $DIRCOUNT "0" InstDirPage
  StrCmp $DIRCOUNT "1" StartDirPage
  InstDirPage:
    StrCpy $DIRCOUNT "1"
    Goto EndDirShow
  StartDirPage:
     StrCpy $DIRCOUNT "2"
     !insertmacro MUI_HEADER_TEXT "Choose WinCoot Start Location" "Choose the folder in which to start WinCoot."
     !insertmacro MUI_INNERDIALOG_TEXT 1041 "WinCoot Start Folder"
     !insertmacro MUI_INNERDIALOG_TEXT 1006 "Setup will make WinCoot icon shortcuts to start in the following folder.$\r$\n$\r$\nTo install in a different folder, click Browse and select another folder. Click Install to continue the installation."
     Goto EndDirShow
  EndDirShow:
FunctionEnd

; for changing wincoot [and gdk-pixbuf-loader] Not sure if here is a good place? Get's a bit messy
Function FinishPagePreFunction
   ; first apply changes to wincoot-for-ccp4i2.bat
   ; this is for ccp4i2 only, not required for normal wincoot.bat
   Var /GLOBAL GUILE_INST_DIR
;   ${WordReplace} "$INSTDIR" "\" "/" "+" $GUILE_INST_DIR
   !insertmacro ReplaceOnLine "yourWinCootdirectory" "$INSTDIR" "5" "$INSTDIR\wincoot-for-ccp4i2.bat.tmp"
;   !insertmacro ReplaceOnLine "yourWinCootdirectoryGuile" "$GUILE_INST_DIR" "6" "$INSTDIR\wincoot.bat.tmp"
   ;we want to change more for Vista (and possibly for Windows 7 too FIXME!)
   ; Maybe not too much any more... (since graphics card issue)
   ${If} ${AtLeastWinVista}
       ; change to run on 1 core only (to enable compositing!)
       !insertmacro AdvReplaceInFile "coot-bin.exe" "start /wait coot-bin.exe" "0" "1" "$INSTDIR\wincoot-for-ccp4i2.bat.tmp"
     ${EndIf}

   ; if we have an old bat file
   ; this we only do for "normal" wincoot.bat; we dont care about i2 for now.
   ${If} $have_bat == "True"
     ; check if wincootbats are different
     Var /Global bat_differ
     StrCpy $bat_differ "False"
     ${TextCompare} "$INSTDIR\wincoot.bat" "$INSTDIR\wincoot.bat.tmp" "FastDiff" "TxtCompResult"

     ${If} $bat_differ == "True"
        ; ask or if silent/update keep old
        Var /GLOBAL keep_old_bat
        StrCpy $keep_old_bat "True"
        ${If} $update = 0
           IfSilent endifbat
            MessageBox MB_ICONQUESTION|MB_YESNO "You already have a (modified) WinCoot batch file (wincoot.bat).$\r$\n\
            Do you want to keep it (dont if you upgrade from <0.9)?" IDYES endifbat
            StrCpy $keep_old_bat "False"
           endifbat:
        ${EndIf}  ; update

        ; replace old with new if requested
        ${If} $keep_old_bat == "False"
             CopyFiles "$INSTDIR\wincoot.bat.tmp" "$INSTDIR\wincoot.bat"
        ${EndIf}

     ${EndIf}  ; bat_differ

  ${Else}
     ; dont have any wincoot.bat
     Rename "$INSTDIR\wincoot.bat.tmp" "$INSTDIR\wincoot.bat"
     Rename "$INSTDIR\wincoot-for-ccp4i2.bat.tmp" "$INSTDIR\wincoot-for-ccp4i2.bat"
  ${EndIf}  ; have_bat nothing further to be one
  ; executable access to everyone
  AccessControl::GrantOnFile /NOINHERIT "$INSTDIR\wincoot.bat" "(BA)" "FullAccess"
  AccessControl::GrantOnFile /NOINHERIT "$INSTDIR\wincoot.bat" "(BU)" "GenericExecute"
  AccessControl::GrantOnFile /NOINHERIT "$INSTDIR\wincoot-for-ccp4i2.bat" "(BA)" "FullAccess"
  AccessControl::GrantOnFile /NOINHERIT "$INSTDIR\wincoot-for-ccp4i2.bat" "(BU)" "GenericExecute"
  AccessControl::GrantOnFile /NOINHERIT "$INSTDIR\bin\dynarama.bat" "(BA)" "FullAccess"
  AccessControl::GrantOnFile /NOINHERIT "$INSTDIR\bin\dynarama.bat" "(BU)" "GenericExecute"

;  for now dont mess with pixbuf query loader
;  ExecWait 'cmd /c ""$INSTDIR\bin\gdk-pixbuf-query-loaders.exe" > "$INSTDIR\etc\gtk-2.0\gdk-pixbuf.loaders""'
  ; set outpath to run where we want to

  SetOutPath "$STARTDIR"
  Call CheckForUpdate
FunctionEnd

Function TxtCompResult
     StrCpy $bat_differ "True"
FunctionEnd

Function MakeIcons
# run only if not update
${If} $update = 0
  SetOutPath "$INSTDIR"
; Shortcuts
  !insertmacro MUI_STARTMENU_WRITE_BEGIN Application
  CreateDirectory "$SMPROGRAMS\$ICONS_GROUP"
  SetOutPath "$STARTDIR"
  CreateShortCut "$QUICKLAUNCH\WinCoot.lnk" "$SYSDIR\cmd.exe" '/c "$INSTDIR\wincoot.bat"' "$INSTDIR\bin\coot-icon.ico" 0 SW_SHOWMINIMIZED
  CreateShortCut "$DESKTOP\WinCoot.lnk" "$SYSDIR\cmd.exe" '/c "$INSTDIR\wincoot.bat"' "$INSTDIR\bin\coot-icon.ico" 0 SW_SHOWMINIMIZED
  CreateShortCut "$SMPROGRAMS\$ICONS_GROUP\WinCoot.lnk" "$SYSDIR\cmd.exe" '/c "$INSTDIR\wincoot.bat"' "$INSTDIR\bin\coot-icon.ico" 0 SW_SHOWMINIMIZED
  Sleep 10
  SetOutPath "$INSTDIR"
  CreateShortCut "$SMPROGRAMS\$ICONS_GROUP\DynaRama.lnk" "$SYSDIR\cmd.exe" '/c "$INSTDIR\bin\dynarama.bat"' "$INSTDIR\bin\rama_all.ico"
  Sleep 10
  WriteIniStr "$INSTDIR\WinCoot.url" "InternetShortcut" "URL" "${PRODUCT_WEB_SITE}"
  WriteIniStr "$INSTDIR\Coot.url" "InternetShortcut" "URL" "${PRODUCT_WEB_SITE_2}"
  CreateShortCut "$SMPROGRAMS\WinCoot\WinCoot Website.lnk" "$INSTDIR\${PRODUCT_NAME}.url"
  Sleep 10
  CreateShortCut "$SMPROGRAMS\WinCoot\Coot Website.lnk" "$INSTDIR\Coot.url"
  Sleep 10
  CreateShortCut "$SMPROGRAMS\WinCoot\Uninstall.lnk" "$INSTDIR\uninst.exe"
  Sleep 10
  CreateShortCut "$SMPROGRAMS\WinCoot\Coot Manual.lnk" "$INSTDIR\doc\coot-user-manual.pdf"
  Sleep 10
  CreateShortCut "$SMPROGRAMS\WinCoot\Coot Keys and Buttons.lnk" "$INSTDIR\doc\crib-sheet.pdf"
  Sleep 10
  CreateShortCut "$SMPROGRAMS\WinCoot\Coot Tutorial.lnk" "$INSTDIR\doc\tutorial.pdf"
  !insertmacro MUI_STARTMENU_WRITE_END
  SetOutPath "$INSTDIR"
 ${Endif}
FunctionEnd

Function ErrorHandler

  Pop $0  ; Abort Flag
  Pop $1  ; Message txt
  Pop $2  ; error code

; convert strings to ints
  IntOp $0 $0 + 0
  IntOp $2 $2 + 0

; the message box could be rather a question to abort!?
  IfSilent +2 0
  MessageBox MB_OK 'Error in Installation. Aborting!$\n$\r$\n$\r$1'

  DetailPrint $1
  SetErrorLevel $2
  ${If} $0 == 1
    Abort $1
  ${EndIf}

FunctionEnd