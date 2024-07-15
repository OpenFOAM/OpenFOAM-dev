/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

#include "bXiProgress.H"
#include "volFields.H"
#include "fvcGrad.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(bXiProgress, 0);
    addToRunTimeSelectionTable(functionObject, bXiProgress, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::bXiProgress::writeFileHeader(const label i)
{
    if (Pstream::master())
    {
        writeHeader(file(), "Combustion progress");
        writeCommented(file(), "Time");

        file() << tab << "progress";

        file() << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::bXiProgress::bXiProgress
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(obr_, name)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::bXiProgress::~bXiProgress()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::bXiProgress::read(const dictionary& dict)
{
    functionObject::read(dict);

    resetName(typeName);

    return true;
}


Foam::wordList Foam::functionObjects::bXiProgress::fields() const
{
    return wordList{"b"};
}


bool Foam::functionObjects::bXiProgress::execute()
{
    return true;
}


bool Foam::functionObjects::bXiProgress::write()
{
    const volScalarField& b =
        mesh_.lookupObject<volScalarField>("b");

    const scalar progress((scalar(1) - b)().weightedAverage(mesh_.V()).value());

    logFiles::write();

    if (Pstream::master())
    {
        writeTime(file());

        file() << tab << progress;

        file() << endl;
    }

    return true;
}


// ************************************************************************* //
