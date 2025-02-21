#!/bin/bash

# This script recursively checks or applies the formatting rules specified in the file ".clang-format". In order to
# check whether all C++ header and source files in the directory "src" are correctly formatted, execute
#     formatting.sh src
# or, equivalently,
#     formatting.sh src	-m check
# In order to format the files, run
#     formatting.sh src -m apply
# You can exclude individual files by adding their names to the file ".clang-format-excludes".
# Note that you need clang in your path for this script to work.
# By providing the -f flag, i.e.
#     formatting.sh src -m check -f
# or
#     formatting.sh src -m apply -f
# the process is done file-wise, i.e. slower but providing more detailed information.

set -o nounset
set -o pipefail

select_clang_format() {
  local version=15
  # Try to find a clang-format binary with version 15
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
  
  if [[ ! -z "${filewise}" ]]; then
    for file in ${files}; do
      unformatted_code=$(echo "${file}" | xargs cat)
      formatted_code=$(echo "${file}" | xargs ${format_binary} -style=file) 
      local replacement_diff=$(diff <(echo "${unformatted_code}") <(echo "${formatted_code}"))

      if [[ ! -z "${replacement_diff}" ]]; then
        echo "Applying formatting to file: " ${file}
        echo "${file}" | xargs ${format_binary} -style=file -i
      fi
    done
  else
    unformatted_code=$(echo "${files}" | xargs cat)
    formatted_code=$(echo "${files}" | xargs ${format_binary} -style=file) 
    local replacement_diff=$(diff <(echo "${unformatted_code}") <(echo "${formatted_code}"))

    if [[ ! -z "${replacement_diff}" ]]; then
      echo "${files}" | xargs ${format_binary} -style=file -i
    fi
  fi

  echo "Finished applying the format specifications."
}

check_formatting () {
  local source_directory="${1}"
  local excludes="${2}"

  local files=$(find_files "${source_directory}" "${excludes}")

  select_clang_format format_binary

  if [[ ! -z "${filewise}" ]]; then
    for file in ${files}; do
      unformatted_code=$(echo "${file}" | xargs cat)
      formatted_code=$(echo "${file}" | xargs ${format_binary} -style=file)
      local replacement_diff=$(diff <(echo "${unformatted_code}") <(echo "${formatted_code}"))

      if [[ ! -z "${replacement_diff}" ]]; then
        echo "ERROR: Source code does not match the format specifications in ${file}"
        exit 1;
      fi
    done
  else
    unformatted_code=$(echo "${files}" | xargs cat)
    formatted_code=$(echo "${files}" | xargs ${format_binary} -style=file)
    local replacement_diff=$(diff <(echo "${unformatted_code}") <(echo "${formatted_code}"))

    if [[ ! -z "${replacement_diff}" ]]; then
      echo "ERROR: Source code does not match the format specifications."
      exit 1;
    fi
  fi
  echo "Finished checking the format specifications."
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

if [[ $# -lt 1 ]]; then
  echo "ERROR: Provide the path to the source files as first arguments!"
  exit 1;
fi
source_directory=$(cd "${1}" && pwd)

# Make sure we are always in the top directory of the repo
relative_script_path=$(dirname "${0}")
absolute_script_path=$(cd ${relative_script_path} && pwd)
cd ${absolute_script_path}/../..

excludes=$(parse_excludes)

mode="check"
filewise=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    -m|--mode)
      mode="$2"
      shift 2
      ;;
    -f|--filewise)
      filewise="true"
      echo "Proceeding filewise"
      shift
      ;;
    *)
      shift
      ;;
  esac
done

case "${mode}" in
  check)
    check_formatting "${source_directory}" "${excludes}"
    ;;
  apply)
    apply_formatting "${source_directory}" "${excludes}"
    ;;
  *)
    echo "Invalid mode. Available options for the mode (-m or --mode) are check and apply."
    exit 1
    ;;
esac

