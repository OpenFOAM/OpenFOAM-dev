/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

//extern "C"
//{
#include "mgridgen.h"
//}

#include "error.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

static const char* notImplementedMessage =
"You are trying to use MGridGen but do not have the MGridGen library loaded.\n"
"This message is from the dummy MGridGen stub library instead.\n"
"\n"
"Normally the MGridGen library will be loaded through the LD_LIBRARY_PATH\n"
"environment variable but you are picking up this dummy library from the\n"
"$FOAM_LIBBIN/dummy directory. Please install MGridGen and make sure the\n"
"libMGridGen.so is in your LD_LIBRARY_PATH.";

#ifdef __cplusplus
extern "C"
#endif
void MGridGen(int, idxtype *, realtype *, realtype *, idxtype *, realtype *,
              int, int, int *, int *, int *, idxtype *)
{
    FatalErrorInFunction
        << notImplementedMessage
        << Foam::exit(Foam::FatalError);
}


// ************************************************************************* //
