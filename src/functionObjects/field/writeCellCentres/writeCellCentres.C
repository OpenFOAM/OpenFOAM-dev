/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2022 OpenFOAM Foundation
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

#include "writeCellCentres.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(writeCellCentres, 0);
    addToRunTimeSelectionTable(functionObject, writeCellCentres, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::writeCellCentres::writeCellCentres
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::writeCellCentres::~writeCellCentres()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::writeCellCentres::execute()
{
    return true;
}


bool Foam::functionObjects::writeCellCentres::write()
{
    volVectorField C
    (
        IOobject
        (
            "C",
            time_.name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_.C(),
        calculatedFvPatchScalarField::typeName
    );

    Log << "    Writing cell-centre field " << C.name()
        << " to " << time_.name() << endl;

    C.write();

    for (direction i=0; i<vector::nComponents; i++)
    {
        volScalarField Ci
        (
            IOobject
            (
                mesh_.C().name() + vector::componentNames[i],
                time_.name(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_.C().component(i)
        );

        Log << "    Writing the "
            << vector::componentNames[i]
            << " component field of the cell-centres " << Ci.name()
            << " to " << time_.name() << endl;

        Ci.write();
    }

    return true;
}


// ************************************************************************* //
