/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2021 OpenFOAM Foundation
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

#include "residuals.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(residuals, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        residuals,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::residuals::residuals
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(obr_, name),
    fieldSet_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::residuals::~residuals()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::residuals::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    dict.lookup("fields") >> fieldSet_;

    resetName(typeName);

    return true;
}


void Foam::functionObjects::residuals::writeFileHeader(const label i)
{
    if (Pstream::master())
    {
        writeHeader(file(), "Residuals");
        writeCommented(file(), "Time");

        forAll(fieldSet_, fieldi)
        {
            const word& fieldName = fieldSet_[fieldi];

            writeFileHeader<scalar>(fieldName);
            writeFileHeader<vector>(fieldName);
            writeFileHeader<sphericalTensor>(fieldName);
            writeFileHeader<symmTensor>(fieldName);
            writeFileHeader<tensor>(fieldName);
        }

        file() << endl;
    }
}


bool Foam::functionObjects::residuals::execute()
{
    return true;
}


bool Foam::functionObjects::residuals::write()
{
    logFiles::write();

    if (Pstream::master())
    {
        writeTime(file());

        forAll(fieldSet_, fieldi)
        {
            const word& fieldName = fieldSet_[fieldi];

            writeResidual<scalar>(fieldName);
            writeResidual<vector>(fieldName);
            writeResidual<sphericalTensor>(fieldName);
            writeResidual<symmTensor>(fieldName);
            writeResidual<tensor>(fieldName);
        }

        file() << endl;
    }

    return true;
}


// ************************************************************************* //
