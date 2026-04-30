using TestItemRunner

# Filter to ensemble-related tests that failed due to SciMLBase v3 prob_func signature change
const FAILING_TESTS = [
    "MultiKetTrajectory construction",
    "MultiKetTrajectory fidelity",
    "SmoothPulseProblem with MultiKetTrajectory",
    "SplinePulseProblem with MultiKetTrajectory",
    "BangBangPulseProblem with MultiKetTrajectory",
    "SplinePulseProblem MultiKet coherent kwarg",
]
@run_package_tests(filter = ti -> ti.name in FAILING_TESTS)
