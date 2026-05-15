using TestItemRunner

# CI matrix shards the testitems across N jobs that run in parallel.
# Locally / when unset, run everything (TOTAL_SHARDS=1).
const TEST_SHARD    = parse(Int, get(ENV, "TEST_SHARD",    "1"))
const TOTAL_SHARDS  = parse(Int, get(ENV, "TOTAL_SHARDS",  "1"))

if TOTAL_SHARDS == 1
    @run_package_tests
else
    # Stable assignment by (filename, name) hash. Same Julia version + same
    # testitem set ⇒ same shard, so coverage is identical across runs in a
    # single CI invocation.
    @run_package_tests filter = ti ->
        mod(hash((ti.filename, ti.name)), TOTAL_SHARDS) == TEST_SHARD - 1
end
