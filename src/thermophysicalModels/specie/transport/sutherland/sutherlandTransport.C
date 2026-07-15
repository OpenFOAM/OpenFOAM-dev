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

#include "sutherlandTransport.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::sutherlandTransport<Thermo>::sutherlandTransport
(
    const Thermo& thermo,
    const dictionary& dict,
    const dictionary& subDict
)
:
    Thermo(thermo),
    As_(subDict.lookup<scalar>
    (
        "As",
        dimensions::dynamicViscosity/sqrt(dimensions::temperature))
    ),
    Ts_(subDict.lookup<scalar>("Ts", dimensions::temperature))
{}


template<class Thermo>
Foam::sutherlandTransport<Thermo>::sutherlandTransport
(
    const Thermo& thermo,
    const dictionary& dict
)
:
    sutherlandTransport(thermo, dict, dict.subDict("transport"))
{}


template<class Thermo>
Foam::sutherlandTransport<Thermo>::sutherlandTransport
(
    const word& name,
    const dictionary& dict
)
:
    sutherlandTransport(Thermo(name, dict), dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::sutherlandTransport<Thermo>::write(Ostream& os) const
{
    Thermo::write(os);

    writeEntry
    (
        os,
        "transport",
        dictionary::entries("As", As_, "Ts", Ts_)
    );
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Thermo>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const sutherlandTransport<Thermo>& st
)
{
    st.write(os);
    return os;
}


// ************************************************************************* //
