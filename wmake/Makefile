#-------------------------------*- makefile -*---------------------------------
# =========                 |
# \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
#  \\    /   O peration     |
#   \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
#    \\/     M anipulation  |
#------------------------------------------------------------------------------
# License
#     This file is part of OpenFOAM.
#
#     OpenFOAM is free software: you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
#     ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#     FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#     for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
#
# File
#     Makefile
#
# Description
#     A generic Makefile, used by wmake
#
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# The Makefile uses a POSIX shell
#------------------------------------------------------------------------------

SHELL      = /bin/sh


#------------------------------------------------------------------------------
# Initialise suffices list
#------------------------------------------------------------------------------

.SUFFIXES:
.SUFFIXES: .o


#------------------------------------------------------------------------------
# Set the directory containing the wmake scripts
#------------------------------------------------------------------------------

WM_SCRIPTS   = $(WM_DIR)/scripts


#------------------------------------------------------------------------------
# declare default paths
#------------------------------------------------------------------------------

LIB_SRC            = $(WM_PROJECT_DIR)/src
LIB_PLATFORMS      = $(WM_PROJECT_DIR)/platforms/$(WM_OPTIONS)/lib
OBJECTS_DIR        = $(MAKE_DIR)/$(WM_OPTIONS)
CLASSES_DIR        = $(MAKE_DIR)/classes

SYS_INC            =
SYS_LIBS           =

PROJECT_INC        = -I$(LIB_SRC)/$(WM_PROJECT)/lnInclude -I$(LIB_SRC)/OSspecific/$(WM_OSTYPE)/lnInclude
PROJECT_LIBS       = -l$(WM_PROJECT)
##
## unused:
## PROJECT_VERSION    = $(shell printf -- "-DPROJECT_VERSION=%d%03d%03d" $$(echo $(WM_PROJECT_VERSION) | sed -e 's/\./ /g'))
##

EXE_INC            =
EXE_LIBS           =

LIB_LIBS           =


#------------------------------------------------------------------------------
# declare all .o files depend on $(OBJECTS_DIR)/options
#------------------------------------------------------------------------------

EXE_DEP         = $(OBJECTS_DIR)/options


#------------------------------------------------------------------------------
# declare default name of libraries and executables
#------------------------------------------------------------------------------

# Library
LIB             = libNULL

# Shared library extension
SO              = so

# Project executable
EXE             = $(WM_PROJECT).out

# Standalone executable
SEXE            = a.out


#------------------------------------------------------------------------------
# set compilation and dependency building rules
#------------------------------------------------------------------------------

GENERAL_RULES = $(WM_DIR)/rules/General
RULES         = $(WM_DIR)/rules/$(WM_ARCH)$(WM_COMPILER)
WMAKE_BIN     = $(WM_DIR)/platforms/$(WM_ARCH)$(WM_COMPILER)

include $(GENERAL_RULES)/general
include $(RULES)/general
include $(RULES)/$(WM_LINK_LANGUAGE)

MAKE_DEP   = @$(MKDEP) $< | $(WM_SCRIPTS)/addCompile $<


#------------------------------------------------------------------------------
# Include PROJECT directory tree file and
# source, object and dependency list files.
# These are constructed by scripts/makeDerivedFiles
#------------------------------------------------------------------------------

include $(OBJECTS_DIR)/options
include $(OBJECTS_DIR)/filesMacros
include $(OBJECTS_DIR)/sourceFiles
include $(OBJECTS_DIR)/objectFiles
include $(OBJECTS_DIR)/localObjectFiles
include $(OBJECTS_DIR)/dependencyFiles

#------------------------------------------------------------------------------
# set header file include paths
#------------------------------------------------------------------------------

LIB_HEADER_DIRS = \
        $(EXE_INC) \
        -IlnInclude \
        -I. \
        $(PROJECT_INC) \
        $(GINC) \
        $(SYS_INC)


#------------------------------------------------------------------------------
# define link statements for executables
#------------------------------------------------------------------------------

$(EXE):  $(OBJECTS)
	     @$(WM_SCRIPTS)/makeTargetDir $(EXE)
	     $(LINKEXE) $(OBJECTS) -L$(LIB_PLATFORMS) \
	     $(EXE_LIBS) $(PROJECT_LIBS) $(SYS_LIBS) $(LINK_LIBS) $(GLIBS) -o $(EXE)

exe:     $(SEXE)
	     @echo \'$(SEXE)\' is up to date.

$(SEXE): $(OBJECTS)
	     @$(WM_SCRIPTS)/makeTargetDir $(SEXE)
	     $(LINKEXE) $(OBJECTS) $(EXE_LIBS) \
         $(SYS_LIBS) $(LINK_LIBS) $(GLIBS) -o $(SEXE)


#------------------------------------------------------------------------------
# define link statements for libraries
#------------------------------------------------------------------------------

objects: $(OBJECTS)

libso:   $(LIB).$(SO)
	     @echo \'$(LIB).$(SO)\' is up to date.

$(LIB).$(SO): $(OBJECTS)
	     @$(WM_SCRIPTS)/makeTargetDir $(LIB)
	     @rm -f so_locations
	     @cd $(OBJECTS_DIR) ; \
	     $(LINKLIBSO) $(LOCAL_OBJECTS) -L$(LIB_PLATFORMS) $(LIB_LIBS) $(GLIB_LIBS) -o $(LIB).$(SO)

lib:     $(LIB).a
	     @echo \'$(LIB).a\' is up to date.

$(LIB).a: $(OBJECTS)
	     @$(WM_SCRIPTS)/makeTargetDir $(LIB)
	     @rm -f $(LIB).a
	     $(AR) $(ARFLAGS) $(LIB).a $(OBJECTS)
	     $(RANLIB) $(LIB).a

libo:    $(LIB).o
	     @echo \'$(LIB).o\' is up to date.

$(LIB).o: $(OBJECTS)
	     @$(WM_SCRIPTS)/makeTargetDir $(LIB)
	     @rm -f $(LIB).o
	     $(LD) -r -o $(LIB).o $(OBJECTS)

jar:     $(LIB).jar
	     @echo \'$(LIB).jar\' is up to date.

$(LIB).jar: $(OBJECTS)
	     @$(WM_SCRIPTS)/makeTargetDir $(LIB)
	     jar cfm $(LIB).jar $(LIB_LIBS) -C $(CLASSES_DIR) .

#------------------------------------------------------------------------------
# Set rule to create the include directory for libraries
#------------------------------------------------------------------------------

lnInclude/uptodate: $(MAKE_DIR)/files $(MAKE_DIR)/options
	     @rm -rf lnInclude ; wmakeLnInclude . ; touch lnInclude/uptodate


#------------------------------------------------------------------------------
# Set depedency rule and include dependency lists
#------------------------------------------------------------------------------

$(OBJECTS_DIR)/dependencies: $(DEPENDENCIES)
	@echo dependencies up to date > $(OBJECTS_DIR)/dependencies

include $(INCLUDE_DEPS)


#------------------------------------------------------------------------------
