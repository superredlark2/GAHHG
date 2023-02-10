.\checkGPU.exe
if %errorlevel%==0 (.\explist.bat) else (goto:exceptionEnd)
exit 0
:exceptionEnd
exit -1