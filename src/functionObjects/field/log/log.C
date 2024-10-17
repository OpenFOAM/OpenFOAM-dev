/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2024 OpenFOAM Foundation
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

#include "log.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(log, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        log,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::log::calc()
{
    if (foundObject<volScalarField>(fieldName_))
    {
        const volScalarField& x = lookupObject<volScalarField>(fieldName_);

        // Cache the current debug setting for dimensionSet
        const bool dimensionSetDebug = dimensionSet::debug;

        // Switch-off dimension checking if requested
        if (!checkDimensions_)
        {
            dimensionSet::debug = 0;
        }

        store
        (
            resultName_,
            clip_ ? Foam::log(max(x, clipValue_)) : Foam::log(x)
        );

        // Reinstate dimension checking
        if (!checkDimensions_)
        {
            dimensionSet::debug = dimensionSetDebug;
        }

        return true;
    }
    else
    {
        cannotFindObject<volScalarField>(fieldName_);

        return false;
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::log::log
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, typeName),
    clip_(false),
    clipValue_(0),
    checkDimensions_(true)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::log::~log()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::log::read(const dictionary& dict)
{
    if (dict.found("clip"))
    {
        clip_ = true;
        dict.lookup("clip") >> clipValue_;
    }

    checkDimensions_ = dict.lookupOrDefault<Switch>("checkDimensions", true);

    return true;
}


// ************************************************************************* //
