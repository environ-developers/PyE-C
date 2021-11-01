import os

from pathlib import Path
from shutil import rmtree, copy2
from subprocess import check_call

from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext

__author__ = "Materialab Research Group"
__contact__ = "oliviero.andreussi@unt.edu"
__version__ = "0.0.1"
__license__ = "GPL"
name = 'pyec'
description = "pyec: Environ Python interface",
long_description = """pyec turns Environ into an engine for embedding or for any other purpose."""


class MakeBuild(build_ext):
    """
    Build instructions for wrapping Environ's Fortran code. Executes targets
    from Makefile.
    """
    def run(self):
        """Executes build."""
        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        """Builds extension."""
        topdir = './'
        self.build_name = self.build_lib + os.sep + name

        try:
            import multiprocessing as mp
            nprocs = mp.cpu_count()
        except ImportError:
            nprocs = 1

        # cleanup
        if os.path.exists(self.build_temp): rmtree(self.build_temp)
        os.makedirs(self.build_temp)

        # remove *.so files
        for file in Path(self.build_temp).glob('*.so'):
            os.remove(file)

        # copy makefiles to temporary build dir
        makefiles = list(Path(topdir + '/install/').glob('*'))
        makefiles += list(Path(topdir + '/src/').glob('**/*.f90'))

        for makefile in makefiles:
            copy2(makefile, self.build_temp)

        # execute makefiles
        check_call(
            ['make', '-f', 'Makefile'] + ['-j', str(nprocs)],
            cwd=self.build_temp,
            env=os.environ.copy(),
        )

        # make build dirs
        for path in (self.build_lib, self.build_name):
            if not os.path.exists(path): os.makedirs(path)

        # copy files to build dir
        files = Path(self.build_temp + os.sep + name + os.sep).glob('*')

        for file in files:
            copy2(file, self.build_name)

        # copy libraries to lib dir
        for files in Path(self.build_temp).glob('*.so'):
            copy2(files, self.build_lib + os.sep)


extensions_pyec = Extension(
    "pyec._pyec",
    sources=[],
)

ext_modules = [
    extensions_pyec,
]

if __name__ == "__main__":

    setup(
        name=name,
        url='',
        description=description,
        version=__version__,
        author=__author__,
        author_email=__contact__,
        license=__license__,
        long_description=long_description,
        python_requires='>=3.6',
        install_requires=['numpy>=1.18.0', 'f90wrap>=0.2.3'],
        extras_require={
            'mpi': [
                'mpi4py>=3.0.0',
            ],
        },
        packages=find_packages('./'),
        ext_modules=ext_modules,
        cmdclass={"build_ext": MakeBuild},
        classifiers=[
            'Development Status :: 1 - Beta',
            'Intended Audience :: Science/Research',
            'Programming Language :: Python :: 3',
            'Topic :: Scientific/Engineering :: Chemistry',
            'Topic :: Scientific/Engineering :: Physics'
        ],
    )
