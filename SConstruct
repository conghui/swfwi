# usage:
# scons [-j n]: compile in release mode
# scons debug=1 [-j n]: compile in debug mode
#
# author: heconghui@gmail.com

import os

# host compiler options
compiler_set = 'gnu'
if compiler_set == 'gnu':
  c_compiler      = ["mpicc",  "-cc=gcc",  "-fopenmp"]
  cxx_compiler    = ["mpicxx", "-cxx=g++", "-fopenmp"]
  linker          = cxx_compiler
  warn_flags      = ["-Wall", "-Wextra", "-Wno-write-strings"]
  optimize_flags  = ["-O2"]
  debug_flags     = ["-O0", "-g"]
  other_flags     = ["-DMPICH_IGNORE_CXX_SEEK", "-DBOOST_LOG_DYN_LINK"]
  link_flags      = ["-O1"]

elif compiler_set == 'intel':
  c_compiler      = ["mpicc",  "-cc=icc",   "-openmp", "-mkl"]
  cxx_compiler    = ["mpicxx", "-cxx=icpc", "-openmp", "-mkl"]
  linker          = cxx_compiler
  warn_flags      = ["-Wall"]
  optimize_flags  = ["-O2"]
  debug_flags     = ["-O0", "-g"]
  other_flags     = ["-DMPICH_IGNORE_CXX_SEEK", "-DBOOST_LOG_DYN_LINK"]
  link_flags      = ["-O1"]

else:
  print "Only GNU or Intel Compiler are supported now"
  Exit(-1)

# set up boost direcotry
boost_root = os.environ['INSTALL_ROOT'] + '/boost/'
boost_inc = boost_root + 'include'
boost_lib = boost_root + 'lib'

mkl_root = os.environ['INSTALL_ROOT'] + '/intel/mkl/'
mkl_inc = mkl_root + 'include'
mkl_lib = mkl_root + 'lib/intel64'

fftw_inc = mkl_inc + '/fftw'

mkl_root = os.environ['INSTALL_ROOT'] + '/softs/install/intel/mkl/'

# set the sub directories (key, value), where value is the name of directory
# please make sure the source code are in src subdirectory
dirlist = [
   ('lib', 'lib'),
   ('bin', 'bin'),
   ('mpi-fwi', 'src/mpi-fwi2d'),
   ('ess-fwi', 'src/ess-fwi2d'),
   ('serial-fwi', 'src/serial-fwi2d'),
   ('fm2d', 'src/fm2d'),
   ('common', 'src/common'),
   ('modeling', 'src/modeling'),
   ('util', 'src/util'),
   ('mdlib', 'src/mdlib'),
]
dirs = dict(dirlist)

# normally, you don't need to modify from the line below
# ---------------------------------------------------------------------- #
is_debug_mode = ARGUMENTS.get('debug', 1)
if int(is_debug_mode):
  print "Debug mode"
  cur_cflags = debug_flags + warn_flags + other_flags
else:
  print "Release mode"
  cur_cflags = optimize_flags + warn_flags + other_flags

# add boost include_dir and lib dir
cur_cflags += ["-isystem", boost_inc, "-isystem", fftw_inc]
libpath     = ["#" + dirs['lib'], boost_lib]
libs        = ["boost_system",    "boost_filesystem",     "boost_thread",
               "boost_date_time", "boost_chrono",         "boost_log_setup",
               "boost_log",       "boost_program_options","boost_timer"]

# add additional includes and libs for gnu compiler
if compiler_set == 'gnu':
  cur_cflags += ["-isystem", mkl_inc]
  libpath    += [mkl_lib]
  libs       += ["mkl_intel_lp64", "mkl_core", "mkl_gnu_thread"]

# setup environment
env = Environment(CC      = c_compiler,
                  CXX     = cxx_compiler,
                  LINK    = linker + link_flags,
                  CCFLAGS = cur_cflags,
                  LIBPATH = libpath,
                  LIBS    = libs,
                  ENV     = os.environ)

# compile
for d in dirlist:
  if not "src/" in d[1]:
    continue

  SConscript(d[1] + "/SConscript",
             variant_dir = d[1].replace('src', 'build'),
             duplicate = 0,
             exports = "env dirs")

