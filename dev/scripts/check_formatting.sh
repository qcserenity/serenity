#!/bin/bash

# This script recursively checks or applies the formatting rules specified in the file ".clang-format". In order to
# check whether all C++ header and source files in the directory "src" are correctly formatted, execute
#     check_formatting.sh src
# In order to format the files, run
#     apply_formatting.sh src
# You can exclude individual files by adding their names to the file ".clang-format-excludes".
# Note that you need clang in your path for this script to work.

set -o errexit
set -o nounset
set -o pipefail

select_clang_format() {
  local version=8
  # Try to find a clang-format binary with version 8
  local binary=$(which clang-format-$version 2> /dev/null)
  if [[ -z ${binary} ]]; then
    binary=$(which clang-format 2> /dev/null)
    if [[ -z ${binary} ]]; then
      echo "Could not find clang-format-${version} or a clang-format binary. Check your PATH"
      exit 1
    fi
  fi

  local tmp=$($binary --version | grep "version $version")
  if [[ -z ${tmp} ]]; then
    echo "Binary '${binary}' does not have required version ${version}"
    exit 1
  fi

  local __resultvar=$1
  eval $__resultvar="'${binary}'"
}

find_files () {
  local source_directory="${1}"
  local excludes="${2}"

  local files=$(find ${source_directory} \( -iname *.h -o -iname *.cpp \) ${excludes})

  echo "${files}"
}

apply_formatting () {
  local source_directory="${1}"
  local excludes="${2}"

  local files=$(find_files "${source_directory}" "${excludes}")

  select_clang_format format_binary

  for file in ${files}; do
    unformatted_code=$(echo "${file}" | xargs cat)
    formatted_code=$(echo "${file}" | xargs ${format_binary} -style=file)
    local replacement_diff=$(diff <(echo "${unformatted_code}") <(echo "${formatted_code}"))

    if [[ ! -z "${replacement_diff}" ]]; then
      echo "Applying formatting to file: " ${file}
      echo "${file}" | xargs ${format_binary} -style=file -i
    fi
  done
}

check_formatting () {
  local source_directory="${1}"
  local excludes="${2}"

  local files=$(find_files "${source_directory}" "${excludes}")

  select_clang_format format_binary
  unformatted_code=$(echo "${files}" | xargs cat)
  formatted_code=$(echo "${files}" | xargs ${format_binary} -style=file)
  local replacement_diff=$(diff <(echo "${unformatted_code}") <(echo "${formatted_code}"))

  if [[ ! -z "${replacement_diff}" ]]; then
    echo "ERROR: Source code does not match the format specifications."
    exit 1;
  fi
}

parse_excludes () {
  local excludes_file=".clang-format-excludes"
  local excludes=""

  if [[ -f ${excludes_file} ]]; then
    local line_count=0
    while read line; do
      line_count=$((line_count + 1))
      if [[ "${line_count}" == "1" ]]; then
        excludes="${excludes} -a -not ( -iname ${line}"
      else
        excludes="${excludes} -or -iname ${line}"
      fi

      if [[ "${line_count} ${excludes_file}" == "$(wc -l ${excludes_file})" ]]; then
        excludes="${excludes} )"
      fi
    done <${excludes_file}
  fi

  echo ${excludes}
}

source_directory=$(cd "${1}" && pwd)

# Make sure we are always in the top directory of the repo
relative_script_path=$(dirname "${0}")
absolute_script_path=$(cd ${relative_script_path} && pwd)
cd ${absolute_script_path}/../..

excludes=$(parse_excludes)

script_name=$(basename "${0}")
if [[ "${script_name}" == "check_formatting.sh" ]]; then
  check_formatting "${source_directory}" "${excludes}"
elif [[ "${script_name}" == "apply_formatting.sh" ]]; then
  apply_formatting "${source_directory}" "${excludes}"
else
  echo "ERROR: This code should never be reached"
  exit 1
fi

