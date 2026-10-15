param(
    [string]$PythonExe,
    [Parameter(Mandatory = $true)][string]$LammpsHome,
    [int]$ExpectedVersion = 20260930
)
$ErrorActionPreference = 'Stop'
$taskProjectRoot = Split-Path -Parent $PSScriptRoot

if (-not $PythonExe) {
    # Prefer the verified interpreter. The pinned NumPy wheel supports <=3.12.
    $taskCandidates = @('Python311', 'Python312', 'Python310')
    foreach ($taskCandidate in $taskCandidates) {
        $taskCandidatePath = Join-Path $env:LOCALAPPDATA "Programs\Python\$taskCandidate\python.exe"
        if (Test-Path -LiteralPath $taskCandidatePath) {
            $PythonExe = $taskCandidatePath
            break
        }
    }
}
if (-not $PythonExe -or -not (Test-Path -LiteralPath $PythonExe)) {
    throw 'Specify -PythonExe with an installed 64-bit Python 3.10-3.12 interpreter; 3.11 is recommended.'
}
& $PythonExe -I -c 'import sys,struct; sys.exit(0 if (3,10) <= sys.version_info < (3,13) and struct.calcsize("P") == 8 else 1)'
if ($LASTEXITCODE -ne 0) { throw 'Pinned dependencies require 64-bit Python 3.10-3.12; 3.11 is recommended.' }

$taskInstallation = (Resolve-Path -LiteralPath $LammpsHome).Path
foreach ($taskComponent in @('Python\lammps\__init__.py', 'bin\liblammps.dll')) {
    if (-not (Test-Path -LiteralPath (Join-Path $taskInstallation $taskComponent))) {
        throw "Missing LAMMPS component: $taskComponent"
    }
}
$taskEnvironment = Join-Path $taskProjectRoot '.venv'
$taskEnvironmentPython = Join-Path $taskEnvironment 'Scripts\python.exe'
if (-not (Test-Path -LiteralPath $taskEnvironmentPython)) {
    & $PythonExe -I -m venv $taskEnvironment
    if ($LASTEXITCODE -ne 0) { throw 'Could not create the project Python environment.' }
}
$taskRequestedPython = & $PythonExe -I -c 'import sys; print(sys.executable)'
if ($LASTEXITCODE -ne 0) { throw 'Could not inspect the requested Python interpreter.' }
$taskActualPython = & $taskEnvironmentPython -I -c 'import sys; print(sys._base_executable)'
if ($LASTEXITCODE -ne 0) { throw 'The existing project Python environment is unusable.' }
if ([System.IO.Path]::GetFullPath($taskRequestedPython.Trim()) -ne [System.IO.Path]::GetFullPath($taskActualPython.Trim())) {
    throw "The existing .venv uses $taskActualPython instead of $taskRequestedPython. Use its original interpreter or recreate the environment before switching."
}
$taskLocal = Join-Path $taskProjectRoot '.local'
New-Item -ItemType Directory -Path $taskLocal -Force | Out-Null
$taskConfiguration = @{ installation = $taskInstallation; expected_version = $ExpectedVersion } | ConvertTo-Json
[System.IO.File]::WriteAllText((Join-Path $taskLocal 'lammps.json'), $taskConfiguration, [System.Text.UTF8Encoding]::new($false))
Write-Output 'Project environment configured. Running the LAMMPS check...'
& (Join-Path $taskProjectRoot 'check-lammps.ps1') -LammpsHome $taskInstallation -ExpectedVersion $ExpectedVersion
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
Write-Output 'Installing the molecular racing dependencies...'
& $taskEnvironmentPython -E -B -m pip install -r (Join-Path $taskProjectRoot 'requirements.txt')
if ($LASTEXITCODE -ne 0) { throw 'Could not install the project dependencies.' }
if (-not (Test-Path -LiteralPath (Join-Path $taskProjectRoot 'frontend\vendor\three.module.js'))) {
    & $taskEnvironmentPython -E -B (Join-Path $taskProjectRoot 'scripts\fetch_frontend.py')
    if ($LASTEXITCODE -ne 0) { throw 'Could not install the browser assets.' }
}
Write-Output 'Ready. Start the prototype with .\start-game.ps1'
exit 0
