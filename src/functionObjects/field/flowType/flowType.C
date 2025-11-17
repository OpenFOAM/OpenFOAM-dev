/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2025 OpenFOAM Foundation
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

#include "flowType.H"
#include "fvcGrad.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(flowType, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        flowType,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::flowType::flowType
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeLocalObjects(obr_, false),
    fieldName_(dict.lookupOrDefault<word>("field", "U")),
    resultName_
    (
        dict.found("result")
      ? dict.lookup<word>("result")
      : name
    )
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::flowType::~flowType()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::flowType::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeLocalObjects::read(dict);

    resetLocalObjectNames(wordList{resultName_, "magSqrDevGradU"});

    return true;
}

Foam::wordList Foam::functionObjects::flowType::fields() const
{
    return wordList(1, fieldName_);
}

bool Foam::functionObjects::flowType::execute()
{
    if (foundObject<volVectorField>(fieldName_))
    {
        const volVectorField& U = lookupObject<volVectorField>(fieldName_);
        const volTensorField devGradU(dev(fvc::grad(U)));

        store
        (
            resultName_,
            -(devGradU && devGradU.T())
            /max
            (
                store("magSqrDevGradU", magSqr(devGradU)),
                dimensionedScalar(sqr(dimRate), rootVSmall)
            )
        );

        return true;
    }
    else
    {
        cannotFindObject<volVectorField>(fieldName_);

        return false;
    }
}

bool Foam::functionObjects::flowType::write()
{
    return writeLocalObjects::write();
}

// ************************************************************************* //
