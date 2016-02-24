#!/bin/bash

# This script is used to format the C/C++ code in a predefine style
# Please run this script before you commit the code.
#
# author: heconghui@gmail.com
#

#--break-blocks \
format_cmd="astyle \
--style=google \
--indent=spaces=2 \
--indent-switches \
--indent-col1-comments \
--indent-preprocessor \
--max-instatement-indent=80 \
--pad-oper \
--pad-header \
--add-brackets \
--convert-tabs \
--align-pointer=name \
--align-reference=name \
--suffix=none \
--preserve-date \
--formatted \
--lineend=linux"

function format_dir()
{
  DIR=$1
  find $DIR -type f \( -name "*.c"    -o -name "*.cpp" -o \
                    -name "*.cxx"  -o -name "*.CC"  -o \
                    -name "*.h"    -o -name "*.hpp" -o \
                    -name "*.hxx"  -o \
                    -name "*.cu"  -o \
                    -name "*.java" \) \
  -print0 | xargs -0 -I {} $format_cmd {}
}

if [[ $# -eq 0 ]]; then
  format_dir .
else
  # iterate all the parameter and format them
  for file in $@; do
    if [[ -f $file ]]; then
      # it is a file
      $format_cmd $file

    elif [[ -d $file ]]; then
      # it is a directory, format all the files in that directory
      format_dir $file
    fi
  done
fi
