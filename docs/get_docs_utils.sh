#!/bin/bash

set -euo pipefail

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$SCRIPT_DIR/.."

WORKFLOW_FILE="$PROJECT_ROOT/.github/workflows/docs.yml"

# Check if workflow file exists
if [[ ! -f "$WORKFLOW_FILE" ]]; then
    echo "GitHub workflow file not found at: $WORKFLOW_FILE"
    exit 1
fi

DOC_TEMPLATE_VERSION=$(grep -E '^\s*DOC_TEMPLATE_VERSION:' "$WORKFLOW_FILE" | sed -E 's/.*DOC_TEMPLATE_VERSION:\s*"([^"]+)".*/\1/')
if [[ -z "$DOC_TEMPLATE_VERSION" ]]; then
    echo "Could not extract DOC_TEMPLATE_VERSION from $WORKFLOW_FILE"
    echo "Expected format: DOC_TEMPLATE_VERSION: \"<version tag here>\""
    exit 1
fi

# Temporary directory for cloning
TEMP_DIR="$PROJECT_ROOT/doc_template_temp"
DOCS_DIR="$PROJECT_ROOT/docs"
UTILS_TARGET="$DOCS_DIR/utils.jl"

# Clean up any existing temporary directory
if [[ -d "$TEMP_DIR" ]]; then
    rm -rf "$TEMP_DIR"
fi

# Clone the repository
echo "Cloning harmoniqs/doc_template at version $DOC_TEMPLATE_VERSION"
# Suppress detached HEAD warning by redirecting stderr and filtering out the warning
if ! git clone --quiet --depth 1 --branch "$DOC_TEMPLATE_VERSION" "https://github.com/harmoniqs/doc_template.git" "$TEMP_DIR" >/dev/null 2>&1; then
    echo "Failed to clone repository at version $DOC_TEMPLATE_VERSION"
    exit 1
fi


UTILS_SOURCE="$TEMP_DIR/utils.jl"

cp "$UTILS_SOURCE" "$UTILS_TARGET"
rm -rf "$TEMP_DIR"

echo "Successfully updated docs/utils.jl with version $DOC_TEMPLATE_VERSION"
