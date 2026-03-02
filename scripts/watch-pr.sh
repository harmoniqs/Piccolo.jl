#!/usr/bin/env bash
# watch-pr.sh — polls a PR's CI checks, re-runs flaky failures, notifies when ready to merge.
#
# Usage: ./scripts/watch-pr.sh [PR_NUMBER] [POLL_INTERVAL_SECONDS] [MAX_RERUNS]
#
# Defaults: PR=97, interval=60s, max_reruns=3

set -euo pipefail

PR="${1:-97}"
INTERVAL="${2:-60}"
MAX_RERUNS="${3:-3}"

declare -A rerun_count  # run_id -> number of times we've triggered a rerun

log() { printf '[%s] %s\n' "$(date '+%H:%M:%S')" "$*"; }

notify() {
    local title="$1" body="$2"
    printf '\a'  # terminal bell
    printf '\n%s\n  %s\n%s\n\n' "$(printf '─%.0s' {1..55})" "$body" "$(printf '─%.0s' {1..55})"
    command -v notify-send &>/dev/null && notify-send "$title" "$body" 2>/dev/null || true
}

branch=$(gh pr view "$PR" --json headRefName --jq '.headRefName')
pr_url=$(gh pr view "$PR" --json url --jq '.url')

log "Watching PR #$PR (${branch})"
log "  URL: ${pr_url}"
log "  Polling every ${INTERVAL}s, max reruns per run: ${MAX_RERUNS}"
log ""

while true; do
    log "Fetching checks..."

    checks=$(gh pr checks "$PR" 2>/dev/null) || {
        log "  gh pr checks failed, retrying in ${INTERVAL}s..."
        sleep "$INTERVAL"
        continue
    }

    if [[ -z "$checks" ]]; then
        log "  No checks found yet. Waiting..."
        sleep "$INTERVAL"
        continue
    fi

    # Columns are tab-separated: name <TAB> status <TAB> duration <TAB> url
    n_pending=$(awk -F'\t' '{print $2}' <<< "$checks" | grep -cE '^(pending|in_progress)$' || true)
    n_fail=$(awk    -F'\t' '{print $2}' <<< "$checks" | grep -c  '^fail$'                  || true)
    n_pass=$(awk    -F'\t' '{print $2}' <<< "$checks" | grep -cE '^(pass|skipping)$'       || true)
    total=$(wc -l <<< "$checks" | tr -d ' ')

    log "  ${n_pass} passed  ${n_fail} failed  ${n_pending} pending  (${total} total)"

    if (( n_pending > 0 )); then
        log "  Checks still running, sleeping ${INTERVAL}s..."
        sleep "$INTERVAL"
        continue
    fi

    if (( n_fail == 0 )); then
        notify "Piccolo CI ✅" "PR #${PR} — all ${total} checks passing. Ready to merge!"
        log "✅ All checks passing — ready to merge!"
        log "   ${pr_url}"
        exit 0
    fi

    # ── There are failures ──────────────────────────────────────────────────
    log "  ${n_fail} check(s) failed. Finding completed failed runs..."

    failed_run_ids=$(gh run list \
        --branch "$branch" \
        --json databaseId,conclusion,status \
        --jq '[.[] | select(.status == "completed" and
               (.conclusion == "failure" or
                .conclusion == "timed_out" or
                .conclusion == "cancelled")) | .databaseId] | .[]' \
        2>/dev/null || true)

    if [[ -z "$failed_run_ids" ]]; then
        log "  No completed failed runs found yet (may still be transitioning). Retrying in ${INTERVAL}s..."
        sleep "$INTERVAL"
        continue
    fi

    gave_up=false
    while IFS= read -r run_id; do
        [[ -z "$run_id" ]] && continue

        count="${rerun_count[$run_id]:-0}"

        if (( count >= MAX_RERUNS )); then
            log "  ❌ Run ${run_id} has failed ${MAX_RERUNS} times — giving up on this run."
            notify "Piccolo CI ❌" "PR #${PR}: run ${run_id} failed ${MAX_RERUNS} times. Manual fix needed."
            gave_up=true
        else
            rerun_count[$run_id]=$(( count + 1 ))
            log "  ↩  Rerunning failed jobs in run ${run_id} (attempt $(( count + 1 ))/${MAX_RERUNS})..."
            if gh run rerun "$run_id" --failed 2>/dev/null; then
                log "     → Rerun triggered."
            else
                log "     → gh run rerun returned non-zero (run may already be rerunning or not rerunnable)."
            fi
        fi
    done <<< "$failed_run_ids"

    if $gave_up; then
        log "One or more runs exceeded the rerun limit. Exiting."
        exit 1
    fi

    log "  Waiting ${INTERVAL}s for reruns to spin up..."
    sleep "$INTERVAL"
done
