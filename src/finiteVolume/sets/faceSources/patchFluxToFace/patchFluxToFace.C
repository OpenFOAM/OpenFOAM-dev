/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2022 OpenFOAM Foundation
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

#include "patchFluxToFace.H"
#include "polyMesh.H"
#include "faceSet.H"
#include "Time.H"
#include "IFstream.H"
#include "fieldDictionary.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(patchFluxToFace, 0);
    addToRunTimeSelectionTable(topoSetSource, patchFluxToFace, word);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::patchFluxToFace::applyToSet
(
    const topoSetSource::setAction action,
    const scalarField& patchFluxField,
    topoSet& set
) const
{
    const polyPatch& patch = mesh().boundaryMesh()[patchName_];

    if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
    {
        Info<< "    Adding all " << fieldName_
            << (inflow_ ? " inflow" : " outflow") << " faces" << endl;

        forAll(patch, facei)
        {
            if
            (
                (inflow_ && patchFluxField[facei] < 0)
             || (!inflow_ && patchFluxField[facei] >= 0)
            )
            {
                set.insert(patch.start() + facei);
            }
        }
    }
    else if (action == topoSetSource::DELETE)
    {
        Info<< "    Removing all " << fieldName_
            << (inflow_ ? " inflow" : " outflow") << " faces" << endl;

        forAll(patch, facei)
        {
            if
            (
                (inflow_ && patchFluxField[facei] < 0)
             || (!inflow_ && patchFluxField[facei] >= 0)
            )
            {
                set.erase(patch.start() + facei);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchFluxToFace::patchFluxToFace
(
    const polyMesh& mesh,
    const word& fieldName,
    const word& patchName,
    const bool inflow
)
:
    topoSetSource(mesh),
    fieldName_(fieldName),
    patchName_(patchName),
    inflow_(inflow)
{}


Foam::patchFluxToFace::patchFluxToFace
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    fieldName_(dict.lookup("field")),
    patchName_(dict.lookup("patch")),
    inflow_(dict.lookup<Switch>("inflow"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::patchFluxToFace::~patchFluxToFace()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::patchFluxToFace::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    typeIOobject<surfaceScalarField> fieldObject
    (
        fieldName_,
        mesh().time().name(),
        mesh(),
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE,
        false
    );

    if (fieldObject.headerOk())
    {
        IFstream str(fieldObject.filePath());

        // Read dictionary
        const fieldDictionary fieldDict
        (
            fieldObject,
            fieldObject.headerClassName()
        );

        const scalarField patchFluxField
        (
            "value",
            fieldDict.subDict("boundaryField").subDict(patchName_),
            mesh().boundaryMesh()[patchName_].size()
        );

        applyToSet(action, patchFluxField, set);
    }
    else
    {
        WarningInFunction
            << "Cannot read flux field " << fieldName_
            << " from time " << mesh().time().name() << endl;
    }
}


// ************************************************************************* //
