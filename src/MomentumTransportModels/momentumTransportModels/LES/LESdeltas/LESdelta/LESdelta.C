/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "LESdelta.H"
#include "calculatedFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(LESdelta, 0);
    defineRunTimeSelectionTable(LESdelta, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LESdelta::LESdelta
(
    const word& name,
    const momentumTransportModel& turbulence
)
:
    momentumTransportModel_(turbulence),
    delta_
    (
        IOobject
        (
            name,
            turbulence.mesh().time().name(),
            turbulence.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        turbulence.mesh(),
        dimensionedScalar(name, dimLength, small),
        calculatedFvPatchScalarField::typeName
    )
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::LESdelta> Foam::LESdelta::New
(
    const word& name,
    const momentumTransportModel& turbulence,
    const dictionary& dict
)
{
    const word deltaType(dict.lookup("delta"));

    Info<< indent
        << "Selecting LES delta type " << deltaType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(deltaType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown LESdelta type "
            << deltaType << nl << nl
            << "Valid LESdelta types are :" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    Info<< incrIndent;
    autoPtr<LESdelta> modelPtr(cstrIter()(name, turbulence, dict));
    Info<< decrIndent;
    return modelPtr;
}


Foam::autoPtr<Foam::LESdelta> Foam::LESdelta::New
(
    const word& name,
    const momentumTransportModel& turbulence,
    const dictionary& dict,
    const dictionaryConstructorTable& additionalConstructors
)
{
    const word deltaType(dict.lookup("delta"));

    Info<< indent
        << "Selecting LES delta type " << deltaType << endl;

    // First on additional ones
    dictionaryConstructorTable::const_iterator cstrIter =
        additionalConstructors.find(deltaType);

    if (cstrIter != additionalConstructors.end())
    {
        Info<< incrIndent;
        autoPtr<LESdelta> modelPtr(cstrIter()(name, turbulence, dict));
        Info<< decrIndent;
        return modelPtr;
    }
    else
    {
        dictionaryConstructorTable::const_iterator cstrIter =
            dictionaryConstructorTablePtr_->find(deltaType);

        if (cstrIter == dictionaryConstructorTablePtr_->end())
        {
            FatalErrorInFunction
                << "Unknown LESdelta type "
                << deltaType << nl << nl
                << "Valid LESdelta types are :" << endl
                << additionalConstructors.sortedToc()
                << " and "
                << dictionaryConstructorTablePtr_->sortedToc()
                << exit(FatalError);
            return autoPtr<LESdelta>();
        }
        else
        {
            Info<< incrIndent;
            autoPtr<LESdelta> modelPtr(cstrIter()(name, turbulence, dict));
            Info<< decrIndent;
            return modelPtr;
        }
    }
}


// ************************************************************************* //
