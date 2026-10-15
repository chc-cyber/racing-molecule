$ErrorActionPreference = 'Stop'
Set-Location (Split-Path $PSScriptRoot -Parent)
$python = Join-Path (Get-Location) '.venv\Scripts\python.exe'
if (!(Test-Path $python)) { throw 'Create .venv and install requirements.txt first; see RELEASING.md.' }
foreach ($tool in @('bootstrap', 'build', 'export_assets', 'export_fourier', 'export_fourier_atlas', 'install_fourier_atlas', 'package_web', 'verify_package')) {
    & $python -E -B "wasm/tools/$tool.py"
    if ($LASTEXITCODE -ne 0) { throw "Web build failed at $tool" }
}
