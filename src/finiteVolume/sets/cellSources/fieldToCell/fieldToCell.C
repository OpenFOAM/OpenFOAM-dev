/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "fieldToCell.H"
#include "polyMesh.H"
#include "cellSet.H"
#include "Time.H"
#include "IFstream.H"
#include "fieldDictionary.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fieldToCell, 0);
    addToRunTimeSelectionTable(topoSetSource, fieldToCell, word);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fieldToCell::applyToSet
(
    const topoSetSource::setAction action,
    const scalarField& field,
    topoSet& set
) const
{
    Info<< "    Field min:" << min(field)
        << " max:" << max(field) << endl;

    if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
    {
        Info<< "    Adding all cells with value of field " << fieldName_
            << " within range " << min_ << ".." << max_ << endl;

        forAll(field, celli)
        {
            if (field[celli] >= min_ && field[celli] <= max_)
            {
                set.insert(celli);
            }
        }
    }
    else if (action == topoSetSource::DELETE)
    {
        Info<< "    Removing all cells with value of field " << fieldName_
            << " within range " << min_ << ".." << max_ << endl;

        forAll(field, celli)
        {
            if (field[celli] >= min_ && field[celli] <= max_)
            {
                set.erase(celli);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fieldToCell::fieldToCell
(
    const polyMesh& mesh,
    const word& fieldName,
    const scalar min,
    const scalar max
)
:
    topoSetSource(mesh),
    fieldName_(fieldName),
    min_(min),
    max_(max)
{}


Foam::fieldToCell::fieldToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    fieldName_(dict.lookup("field")),
    min_(dict.lookup<scalar>("min")),
    max_(dict.lookup<scalar>("max"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fieldToCell::~fieldToCell()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fieldToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    typeIOobject<volScalarField> fieldObject
    (
        fieldName_,
        mesh().time().name(),
        mesh(),
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE,
        false
    );

    if (!fieldObject.IOobject::headerOk())
    {
        WarningInFunction
            << "Cannot read field " << fieldName_
            << " from time " << mesh().time().name() << endl;
    }
    else if (fieldObject.IOobject::headerOk())
    {
        IFstream str(fieldObject.filePath());

        // Read dictionary
        const fieldDictionary fieldDict
        (
            fieldObject,
            volScalarField::typeName
        );

        const scalarField internalVals
        (
            "internalField",
            fieldDict, mesh().nCells()
        );

        applyToSet(action, internalVals, set);
    }
    else
    {
        typeIOobject<volVectorField> fieldObject
        (
            fieldName_,
            mesh().time().name(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE,
            false
        );

        if (fieldObject.IOobject::headerOk())
        {
            IFstream str(fieldObject.filePath());

            // Read dictionary
            const fieldDictionary fieldDict
            (
                fieldObject,
                volVectorField::typeName
            );

            const vectorField internalVals
            (
                "internalField",
                fieldDict, mesh().nCells()
            );

            applyToSet(action, mag(internalVals), set);
        }
        else
        {
            WarningInFunction
            << "Cannot handle fields of type " << fieldObject.headerClassName()
            << endl;
        }
    }
}


// ************************************************************************* //
