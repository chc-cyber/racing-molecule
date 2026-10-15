param([int]$Port = 8765)
$ErrorActionPreference = 'Stop'
Invoke-RestMethod -Uri "http://127.0.0.1:$Port/api/shutdown" -Method Post -ContentType 'application/json' -Body '{}' | Out-Null
Write-Output 'Shutdown requested for the local molecular racing server.'
