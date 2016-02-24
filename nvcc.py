"""SCons.Tool.nvcc

Tool-specific initialization for NVIDIA CUDA Compiler.

There normally shouldn't be any need to import this module directly.
It will usually be imported through the generic SCons.Tool.Tool()
selection method.

"""

import SCons.Tool
import SCons.Scanner.C
import SCons.Defaults

CUDASuffixes = ['.cu']

# make a CUDAScanner for finding #includes
# cuda uses the c preprocessor, so we can use the CScanner
CUDAScanner = SCons.Scanner.C.CScanner()

def which(program):
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

def add_common_nvcc_variables(env):
  """
  Add underlying common "NVIDIA CUDA compiler" variables that
  are used by multiple builders.
  """

  # "NVCC common command line"
  if not env.has_key('_NVCCCOMCOM'):
    # prepend -Xcompiler before each flag

    # these flags are common to both static and shared compilations
    env['_NVCCCOMCOM'] = '${_concat("-Xcompiler ", CPPFLAGS, "", __env__)} $_CPPDEFFLAGS $_CPPINCFLAGS'

    # wrap up all these environment variables inside -Xcompiler ""
    env['_NVCCWRAPCFLAGS'] =     '${_concat("-Xcompiler ", CFLAGS,     "", __env__)}'
    env['_NVCCWRAPSHCFLAGS'] =   '${_concat("-Xcompiler ", SHCFLAGS,   "", __env__)}'
    env['_NVCCWRAPCCFLAGS'] =    '${_concat("-Xcompiler ", CCFLAGS,    "", __env__)}'
    env['_NVCCWRAPSHCCFLAGS'] =  '${_concat("-Xcompiler ", SHCCFLAGS,  "", __env__)}'
    # XXX should these be wrapped as well?  not sure -jph
    #env['_NVCCWRAPCXXFLAGS'] =   '${_concat("-Xcompiler ", CXXFLAGS,   "", __env__)}'
    #env['_NVCCWRAPSHCXXFLAGS'] = '${_concat("-Xcompiler ", SHCXXFLAGS, "", __env__)}'

def generate(env):
  """
  Add Builders and construction variables for CUDA compilers to an Environment.
  """

  static_obj, shared_obj = SCons.Tool.createObjBuilders(env)

  for suffix in CUDASuffixes:
    # Add this suffix to the list of things buildable by Object
    static_obj.add_action('$CUDAFILESUFFIX', '$NVCCCOM')
    shared_obj.add_action('$CUDAFILESUFFIX', '$SHNVCCCOM')
    static_obj.add_emitter(suffix, SCons.Defaults.StaticObjectEmitter)
    shared_obj.add_emitter(suffix, SCons.Defaults.SharedObjectEmitter)

    # Add this suffix to the list of things scannable
    SCons.Tool.SourceFileScanner.add_scanner(suffix, CUDAScanner)

  add_common_nvcc_variables(env)

  # set the "CUDA Compiler Command" environment variable
  env['NVCC'] = 'nvcc'
  env['SHNVCC'] = 'nvcc'

  # set the include path, and pass both c compiler flags and c++ compiler flags
  env['NVCCFLAGS'] = SCons.Util.CLVar('')
  env['SHNVCCFLAGS'] = SCons.Util.CLVar('') + ' -shared'

  # 'NVCC Command'
  env['NVCCCOM']   = '$NVCC -o $TARGET -c $NVCCFLAGS $_NVCCWRAPCFLAGS $_NVCCWRAPCCFLAGS $_NVCCCOMCOM $SOURCES'
  env['SHNVCCCOM'] = '$SHNVCC -o $TARGET -c $SHNVCCFLAGS $_NVCCWRAPSHCFLAGS $_NVCCWRAPSHCCFLAGS $_NVCCCOMCOM $SOURCES'

  # the suffix of CUDA source files is '.cu'
  env['CUDAFILESUFFIX'] = '.cu'

  command  = 'nvcc'
  nvccPath = which(command)
  cudaRoot = nvccPath[:-len(command + 'bin/')];
  exe_path = cudaRoot + 'bin'
  lib_path = cudaRoot + 'lib64'

  env.PrependENVPath('PATH', exe_path)

def exists(env):
  return env.Detect('nvcc')
