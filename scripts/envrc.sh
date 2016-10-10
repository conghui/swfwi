#!/usr/bin/env zsh
# vim: ft=sh
# vim: fdm=marker fdl=0
export VISUAL=${HOME}/.linuxbrew/bin/nvim
function remove_PATH() { # {{{
  # only one path is allowed at a time
  if [[ -d "$1" ]]; then
    export PATH=`echo -n $PATH | awk -v RS=: -v ORS=: '$0 != "'$1'"' | sed 's/:$//'`
  fi
}
# }}}
function prepend_PATH() { # {{{
  # only one path is allowed at a time
  remove_PATH $1
  if [[ -d "$1" ]]; then
    export PATH="$1:$PATH"
  fi
}
# }}}

source /opt/intel/composer_xe_2013_sp1/bin/compilervars.sh intel64
prepend_PATH ~/softs/install/mpich-3.1.3-gnu/bin
remove_PATH ~/.linuxbrew/bin
remove_PATH ~/.linuxbrew/sbin

