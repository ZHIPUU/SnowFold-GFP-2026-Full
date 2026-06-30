$log = 'D:\workspace\round28_local\watcher.log'
New-Item -ItemType Directory -Force 'D:\workspace\round28_local' | Out-Null
"[$(Get-Date)] R28 watcher started; waiting for R27 to finish" | Add-Content $log
while ($true) {
  $p = Get-CimInstance Win32_Process | Where-Object { $_.CommandLine -like '*r27_diverge.py*' }
  if (-not $p) { break }
  "[$(Get-Date)] R27 still running; sleep 5 min" | Add-Content $log
  Start-Sleep -Seconds 300
}
"[$(Get-Date)] R27 finished; launching R28" | Add-Content $log
Set-Location 'D:\workspace'
python -u r28_local_mutscan.py *> 'D:\workspace\round28_local\r28.log'
"[$(Get-Date)] R28 finished" | Add-Content $log
