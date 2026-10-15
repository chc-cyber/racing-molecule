param(
    [string]$LammpsHome,
    [int]$ExpectedVersion
)
$ErrorActionPreference = 'Stop'
$taskPython = Join-Path $PSScriptRoot '.venv\Scripts\python.exe'
if (-not (Test-Path -LiteralPath $taskPython)) {
    throw 'Run scripts/setup.ps1 first to create the project Python environment.'
}
$taskArguments = @('-E', '-B', (Join-Path $PSScriptRoot 'scripts\check_lammps.py'))
if ($LammpsHome) { $taskArguments += @('--lammps-home', $LammpsHome) }
if ($ExpectedVersion) { $taskArguments += @('--expected-version', $ExpectedVersion) }
& $taskPython @taskArguments
exit $LASTEXITCODE
