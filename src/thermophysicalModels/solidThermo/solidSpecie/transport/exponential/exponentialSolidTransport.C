/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

#include "exponentialSolidTransport.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::exponentialSolidTransport<Thermo>::exponentialSolidTransport
(
    const word& name,
    const dictionary& dict,
    const dictionary& subDict
)
:
    Thermo(name, dict),
    kappa0_(subDict.lookup<scalar>("kappa0", dimThermalConductivity)),
    n0_(subDict.lookup<scalar>("n0", dimless)),
    Tref_(subDict.lookup<scalar>("Tref", dimTemperature))
{}


template<class Thermo>
Foam::exponentialSolidTransport<Thermo>::exponentialSolidTransport
(
    const word& name,
    const dictionary& dict
)
:
    exponentialSolidTransport(name, dict, dict.subDict("transport"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::exponentialSolidTransport<Thermo>::exponentialSolidTransport::write
(
    Ostream& os
) const
{
    Thermo::write(os);

    writeEntry
    (
        os,
        "transport",
        dictionary::entries("kappa0", kappa0_, "n0", n0_, "Tref", Tref_)
    );
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Thermo>
Foam::Ostream& Foam::operator<<
(
    Ostream& os, const exponentialSolidTransport<Thermo>& et
)
{
    et.write(os);
    return os;
}


// ************************************************************************* //
