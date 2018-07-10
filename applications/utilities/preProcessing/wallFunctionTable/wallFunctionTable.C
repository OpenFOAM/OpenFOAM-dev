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

Application
    wallFunctionTable

Description
    Generates a table suitable for use by tabulated wall functions.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "tabulatedWallFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    IOdictionary dict
    (
        IOobject
        (
            "wallFunctionDict",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED
        )
    );

    autoPtr<tabulatedWallFunctions::tabulatedWallFunction>
        twf(tabulatedWallFunctions::tabulatedWallFunction::New(dict, mesh));

//    twf->writeData(Info);

    twf->write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
