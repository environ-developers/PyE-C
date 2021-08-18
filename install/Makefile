#----------------------------------------------------------------------------------------
#
# Copyright (C) 2021 ENVIRON (www.quantum-environ.org)
#
#----------------------------------------------------------------------------------------
#
#     This file is part of Environ version 2.0
#     
#     Environ 2.0 is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 2 of the License, or
#     (at your option) any later version.
#     
#     Environ 2.0 is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more detail, either the file
#     `License' in the root directory of the present distribution, or
#     online at <http://www.gnu.org/licenses/>.
#
#----------------------------------------------------------------------------------------
#
# Authors: Edan Bainglass (Department of Physics, UNT)
#
#----------------------------------------------------------------------------------------
#
# Makefile for PyE-C
#
#----------------------------------------------------------------------------------------

ENVDIR := $(or ${envdir}, ../../../)
include ${ENVDIR}/make.inc

MODFLAGS = $(BASEMOD_FLAGS)
ENVLIBS = -L${ENVDIR}/libs/ -lenv -L${ENVDIR}/libs/ -lfft -L${ENVDIR}/libs/ -lutil

PYTHON = python3

ifeq ($(MPIF90), mpiifort)
	FCOMP = intelem
	LIBS = ${BLAS_LIBS}
	F2PYFLAGS = --f90exec=mpiifort
endif

UNAME = $(shell uname)

ifeq (${UNAME}, Darwin)
    LIBTOOL = libtool -static -o
else
    LIBTOOL = ar src
endif

# PROJECT CONFIG
PYTHON_MODN = pyec
KIND_MAP = kind_map

# List all source files required for the project

LIBSRC_SOURCES = setup_interface calc_interface control_interface output_interface

# File names
LIBSRC_FILES = $(addsuffix .f90, ${LIBSRC_SOURCES})

# Object files
LIBSRC_OBJECTS = $(addsuffix .o, ${LIBSRC_SOURCES})

## List all source files that require a Python interface
# names without suffix, f90 sources
LIBSRC_WRAP_SOURCES = setup_interface calc_interface control_interface output_interface

# file names
LIBSRC_WRAP_FILES = $(addsuffix .f90, ${LIBSRC_WRAP_SOURCES})

# object files
LIBSRC_WRAP_OBJECTS = $(addsuffix .o, ${LIBSRC_WRAP_SOURCE})

# Relevant suffixes
.SUFFIXES: .f90 .fpp

.PHONY: all clean

all: _${PYTHON_MODN}.so

clean:
	@ rm -fr *.mod *.a *.so *.fpp f90wrap*.f90 f90wrap*.o *.o src.linux* .f2py_f2cmap .libs pyec/

pyec_interface.o : ${ENVDIR}/utils/base_io.o
pyec_interface.o : ${ENVDIR}/utils/char_ops.o
pyec_interface.o : ${ENVDIR}/utils/io.o

libsrc.a: ${LIBSRC_OBJECTS}
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@
	
_${PYTHON_MODN}.so: libsrc.a ${LIBSRC_FILES}
	f90wrap -m ${PYTHON_MODN} ${LIBSRC_WRAP_FILES} -k ${KIND_MAP} -v -P
	f2py-f90wrap --fcompiler=$(FCOMP) ${LIBS} ${F2PYFLAGS} --build-dir . -c -m _${PYTHON_MODN} -L. -lsrc ${ENVLIBS} f90wrap*.f90

include make.depend
