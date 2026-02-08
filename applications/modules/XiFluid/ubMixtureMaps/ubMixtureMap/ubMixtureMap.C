/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "ubMixtureMap.H"
#include "uMixture.H"
#include "bMixture.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ubMixtureMap, 0);
    defineRunTimeSelectionTable(ubMixtureMap, thermo);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ubMixtureMap::ubMixtureMap
(
    const uRhoMulticomponentThermo& uThermo,
    const bRhoMulticomponentThermo& bThermo
)
:
    uThermo_(uThermo),
    bThermo_(bThermo)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ubMixtureMap::~ubMixtureMap()
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::ubMixtureMap> Foam::ubMixtureMap::New
(
    const uRhoMulticomponentThermo& uThermo,
    const bRhoMulticomponentThermo& bThermo
)
{
    const word mapType
    (
        dynamicCast<const uMixture&>(uThermo).mixtureType()
      + dynamicCast<const bMixture&>(bThermo).mixtureType()
    );

    thermoConstructorTable::iterator cstrIter =
        thermoConstructorTablePtr_->find(mapType);

    if (cstrIter == thermoConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown ubMixtureMap " << mapType << nl << nl
            << "Valid ubMixtureMaps are : " << endl
            << thermoConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<ubMixtureMap>(cstrIter()(uThermo, bThermo));
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// ************************************************************************* //
