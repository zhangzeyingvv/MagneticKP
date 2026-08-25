#!/bin/sh

# Canonicalize and validate every generated authoring source. The .nb files
# under MagneticKP/Documentation are generated artifacts; never edit them.

set -eu

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
project_root=$(dirname "$(dirname "$script_dir")")
codex_root=${CODEX_HOME:-"$HOME/.codex"}
guard=${MATHEMATICA_NOTEBOOK_GUARD:-"$codex_root/skills/mathematica-notebook/scripts/notebook_guard.py"}

if [ ! -f "$guard" ]; then
  echo "Notebook guard not found: $guard" >&2
  echo "Set MATHEMATICA_NOTEBOOK_GUARD to notebook_guard.py." >&2
  exit 2
fi

build_language() {
  source_dir=$1
  expected_dir=$2
  target_dir=$3

  find "$source_dir" -type f -name '*.source.wl' -print | sort |
    while IFS= read -r source_path; do
      relative_path=${source_path#"$source_dir"/}
      relative_base=${relative_path%.source.wl}
      notebook_path="$target_dir/$relative_base.nb"
      expected_path="$expected_dir/$relative_base.expected.txt"

      mkdir -p "$(dirname "$notebook_path")"
      echo "Generating $notebook_path"
      python3 "$guard" canonicalize "$source_path" "$notebook_path"
      python3 "$guard" validate "$notebook_path" --expect-file "$expected_path"
    done
}

build_language \
  "$script_dir/Sources" \
  "$script_dir/Expected" \
  "$project_root/MagneticKP/Documentation/English"

build_language \
  "$script_dir/SourcesChineseSimplified" \
  "$script_dir/ExpectedChineseSimplified" \
  "$project_root/MagneticKP/Documentation/ChineseSimplified"

echo "Generated and validated 28 bilingual MagneticKP notebooks."
