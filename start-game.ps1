param(
    [int]$Port = 8765,
    [string]$Molecule
)
$ErrorActionPreference = 'Stop'
$projectDirectory = $PSScriptRoot
$python = Join-Path $projectDirectory '.venv\Scripts\python.exe'
if (-not (Test-Path -LiteralPath $python)) {
    throw 'Initialize the project with scripts/setup.ps1 first.'
}
Push-Location -LiteralPath $projectDirectory
try {
    $arguments = @('-E', '-B', '-m', 'molecular_car.server', '--port', $Port)
    if ($Molecule) { $arguments += @('--molecule', $Molecule) }
    & $python @arguments
    exit $LASTEXITCODE
} finally {
    Pop-Location
}
