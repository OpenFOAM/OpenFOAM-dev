/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "reactingMultiphaseParcelInjectionData.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(reactingMultiphaseParcelInjectionData, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reactingMultiphaseParcelInjectionData::
reactingMultiphaseParcelInjectionData()
:
    reactingParcelInjectionData(),
    YGas_(),
    YLiquid_(),
    YSolid_()
{}


Foam::reactingMultiphaseParcelInjectionData::
reactingMultiphaseParcelInjectionData
(
    const dictionary& dict
)
:
    reactingParcelInjectionData(dict),
    YGas_(dict.lookup("YGas")),
    YLiquid_(dict.lookup("YLiquid")),
    YSolid_(dict.lookup("YSolid"))
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::reactingMultiphaseParcelInjectionData::
~reactingMultiphaseParcelInjectionData()
{}


// ************************************************************************* //
