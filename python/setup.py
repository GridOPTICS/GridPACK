# -------------------------------------------------------------
# file: setup.py
# -------------------------------------------------------------
# Package metadata lives in pyproject.toml.  This file exists only
# to drive the CMake build of the pybind11 extension
# `gridpack._gridpack`; the rest of the packaging (name, version,
# dependencies, entry points, ...) comes from pyproject.toml under
# the setuptools build backend.
# -------------------------------------------------------------

import os
import sys
import subprocess

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        try:
            subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError(
                "CMake must be installed to build the following extensions: "
                + ", ".join(e.name for e in self.extensions)
            )

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        cmake_args = [
            '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
            '-DPython_EXECUTABLE=' + sys.executable,
            '-DCMAKE_VERBOSE_MAKEFILE=YES',
        ]

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]

        cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
        build_args += ['--', '-j2']

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(
            env.get('CXXFLAGS', ''),
            self.distribution.get_version(),
        )
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(
            ['cmake', ext.sourcedir] + cmake_args,
            cwd=self.build_temp,
            env=env,
        )
        subprocess.check_call(
            ['cmake', '--build', '.'] + build_args,
            cwd=self.build_temp,
        )


setup(
    ext_modules=[CMakeExtension('gridpack._gridpack')],
    cmdclass=dict(build_ext=CMakeBuild),
    scripts=[
        'src/hello.py',
        'src/hello_comm.py',
        'src/task_manager.py',
        'src/hadrec.py',
        'src/hadrec_comm.py',
        'src/dsf.py',
        'src/dsf2.py',
        'src/emt.py',
        'src/pf.py',
        'src/stes.py',
    ],
)
