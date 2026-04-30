using TestItemRunner

# Filter to ensemble-related tests that failed due to SciMLBase v3 breaking changes
const FAILING_TESTS = [
    # MultiKetTrajectory (EnsembleProblem prob_func + EnsembleSolution indexing)
    "MultiKetTrajectory construction",
    "MultiKetTrajectory fidelity",
    "SmoothPulseProblem with MultiKetTrajectory",
    "SplinePulseProblem with MultiKetTrajectory",
    "BangBangPulseProblem with MultiKetTrajectory",
    "SplinePulseProblem MultiKet coherent kwarg",
    # MultiDensityTrajectory (same EnsembleSolution indexing pattern)
    "MultiDensityTrajectory construction",
    "MultiDensityTrajectory callable and indexing",
    "MultiDensityTrajectory fidelity",
    # Rollout paths that use EnsembleProblem / EnsembleSolution
    "rollout - MultiKetTrajectory",
    "Two ways to check fidelity",
]
@run_package_tests(filter = ti -> ti.name in FAILING_TESTS)
