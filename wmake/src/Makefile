#------------------------------------------------------------------------------
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
#     Makefile used to compile wmake utilities.
#
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# The Makefile use a POSIX shell
#------------------------------------------------------------------------------

SHELL      = /bin/sh

#------------------------------------------------------------------------------
# Initialise suffices list
#------------------------------------------------------------------------------

.SUFFIXES:
.SUFFIXES: .o


#------------------------------------------------------------------------------
# set compilation and dependency building rules
#------------------------------------------------------------------------------

GENERAL_RULES = $(WM_DIR)/rules/General
RULES         = $(WM_DIR)/rules/$(WM_ARCH)$(WM_COMPILER)
WMAKE_BIN     = $(WM_DIR)/platforms/$(WM_ARCH)$(WM_COMPILER)

include $(RULES)/general
include $(RULES)/$(WM_LINK_LANGUAGE)


#------------------------------------------------------------------------------
# targets
#------------------------------------------------------------------------------

all:  $(WMAKE_BIN)/dirToString $(WMAKE_BIN)/wmkdep $(WMAKE_BIN)/wmkdepend

clean:
	rm -f $(WMAKE_BIN)/* 2>/dev/null


$(WMAKE_BIN)/dirToString: dirToString.c
	@mkdir -p $(WMAKE_BIN)
	$(cc) $(cFLAGS) dirToString.c -o $(WMAKE_BIN)/dirToString


$(WMAKE_BIN)/wmkdep: wmkdep.l
	@mkdir -p $(WMAKE_BIN)
	flex wmkdep.l
	$(cc) $(cFLAGS) lex.yy.c -o $(WMAKE_BIN)/wmkdep
	@rm -f lex.yy.c 2>/dev/null


# for bootstrapping - use generated files directly (instead of from .atg file)
$(WMAKE_BIN)/wmkdepend: wmkdepend.cpp \
    wmkdependParser.cpp wmkdependScanner.cpp \
    wmkdependParser.h   wmkdependScanner.h
	@mkdir -p $(WMAKE_BIN)
	$(CC) $(c++FLAGS) \
	     wmkdepend.cpp wmkdependParser.cpp wmkdependScanner.cpp \
	    -o $(WMAKE_BIN)/wmkdepend

#------------------------------------------------------------------------------
