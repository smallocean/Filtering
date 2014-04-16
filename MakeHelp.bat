@echo off
REM -- First make map file from Microsoft Visual C++ generated resource.h
echo // MAKEHELP.BAT generated Help Map file.  Used by MCUBE.HPJ. >"hlp\Mcube.hm"
echo. >>"hlp\Mcube.hm"
echo // Commands (ID_* and IDM_*) >>"hlp\Mcube.hm"
makehm ID_,HID_,0x10000 IDM_,HIDM_,0x10000 resource.h >>"hlp\Mcube.hm"
echo. >>"hlp\Mcube.hm"
echo // Prompts (IDP_*) >>"hlp\Mcube.hm"
makehm IDP_,HIDP_,0x30000 resource.h >>"hlp\Mcube.hm"
echo. >>"hlp\Mcube.hm"
echo // Resources (IDR_*) >>"hlp\Mcube.hm"
makehm IDR_,HIDR_,0x20000 resource.h >>"hlp\Mcube.hm"
echo. >>"hlp\Mcube.hm"
echo // Dialogs (IDD_*) >>"hlp\Mcube.hm"
makehm IDD_,HIDD_,0x20000 resource.h >>"hlp\Mcube.hm"
echo. >>"hlp\Mcube.hm"
echo // Frame Controls (IDW_*) >>"hlp\Mcube.hm"
makehm IDW_,HIDW_,0x50000 resource.h >>"hlp\Mcube.hm"
REM -- Make help for Project MCUBE


echo Building Win32 Help files
start /wait hcw /C /E /M "hlp\Mcube.hpj"
if errorlevel 1 goto :Error
if not exist "hlp\Mcube.hlp" goto :Error
if not exist "hlp\Mcube.cnt" goto :Error
echo.
if exist Debug\nul copy "hlp\Mcube.hlp" Debug
if exist Debug\nul copy "hlp\Mcube.cnt" Debug
if exist Release\nul copy "hlp\Mcube.hlp" Release
if exist Release\nul copy "hlp\Mcube.cnt" Release
echo.
goto :done

:Error
echo hlp\Mcube.hpj(1) : error: Problem encountered creating help file

:done
echo.
