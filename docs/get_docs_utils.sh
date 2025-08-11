#!/bin/bash

set -euo pipefail

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$SCRIPT_DIR/.."

# if argument is provided, use it as the DOC_TEMPLATE_VERSION
if [[ $# -gt 0 ]]; then
    DOC_TEMPLATE_VERSION="$1"
else
    WORKFLOW_FILE="$PROJECT_ROOT/.github/workflows/docs.yml"

    # Check if workflow file exists
    if [[ ! -f "$WORKFLOW_FILE" ]]; then
        echo "GitHub workflow file not found at: $WORKFLOW_FILE"
        exit 1
    fi

    DOC_TEMPLATE_VERSION=$(grep -E '^\s*DOC_TEMPLATE_VERSION:' "$WORKFLOW_FILE" | sed -E 's/.*DOC_TEMPLATE_VERSION:\s*"([^"]+)".*/\1/')
fi

if [[ -z "$DOC_TEMPLATE_VERSION" ]]; then
    echo "DOC_TEMPLATE_VERSION is not set"
    echo "Please provide a version tag as an arg or ensure it is set in $WORKFLOW_FILE"
    echo "Could not extract DOC_TEMPLATE_VERSION from $WORKFLOW_FILE"
    echo "Expected format: DOC_TEMPLATE_VERSION: \"<version tag here>\""
    exit 1
fi

# Temporary directory for cloning
TEMP_DIR="$PROJECT_ROOT/doc_template_temp"
DOCS_DIR="$PROJECT_ROOT/docs"

# Clone the repository
echo "Grabbing PiccoloDocsTemplate at version $DOC_TEMPLATE_VERSION"
julia --project="$PROJECT_ROOT/docs" -e "
using Pkg; Pkg.add(url=\"https://github.com/harmoniqs/PiccoloDocsTemplate.jl\", rev=\"$DOC_TEMPLATE_VERSION\")
"

echo "Successfully updated PiccoloDocsTemplate with version $DOC_TEMPLATE_VERSION"
