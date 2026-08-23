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

#include "carrierThermo.H"
#include "physicalProperties.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class CarrierThermo, class Thermo>
Foam::autoPtr<CarrierThermo> Foam::carrierThermo::New
(
    const cloud& c,
    const clouds::carried& carriedCloud,
    const word& phaseName
)
{
    const word thermoName =
        IOobject::groupName(physicalProperties::typeName, phaseName);

    if (!c.mesh().poly().foundObject<Thermo>(thermoName))
    {
        return autoPtr<CarrierThermo>(nullptr);
    }

    const Thermo& thermo =
        c.mesh().poly().lookupObject<Thermo>(thermoName);

    const word thermoTypeName =
        CarrierThermo::derivedThermoName() + "<" + thermo.mixtureName() + ">";

    typename CarrierThermo::cloudConstructorTable::iterator cstrIter =
        CarrierThermo::cloudConstructorTablePtr_->find(thermoTypeName);

    if (cstrIter == CarrierThermo::cloudConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown " << CarrierThermo::typeName << " type " << nl << nl
            << thermoTypeName << nl << nl
            << "Valid " << CarrierThermo::typeName << " types are:" << nl
            << CarrierThermo::cloudConstructorTablePtr_->sortedToc() << nl
            << exit(FatalError);
    }

    return autoPtr<CarrierThermo>(cstrIter()(c, carriedCloud, thermo));
}


// ************************************************************************* //
