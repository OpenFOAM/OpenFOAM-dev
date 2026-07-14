/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2026 OpenFOAM Foundation
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

#include "AndradeTransport.H"
#include "dictionary.H"
#include "units.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::AndradeTransport<Thermo>::AndradeTransport
(
    const word& name,
    const dictionary& dict,
    const dictionary& subDict
)
:
    Thermo(name, dict),
    muCoeffs_(subDict.lookup<coeffList>("muCoeffs", units::none)),
    kappaCoeffs_(subDict.lookup<coeffList>("kappaCoeffs", units::none))
{}


template<class Thermo>
Foam::AndradeTransport<Thermo>::AndradeTransport
(
    const word& name,
    const dictionary& dict
)
:
    AndradeTransport(name, dict, dict.subDict("transport"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::AndradeTransport<Thermo>::write(Ostream& os) const
{
    Thermo::write(os);

    writeEntry
    (
        os,
        "transport",
        dictionary::entries
        (
            "muCoeffs", muCoeffs_,
            "kappaCoeffs", kappaCoeffs_
        )
    );
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Thermo>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const AndradeTransport<Thermo>& pt
)
{
    pt.write(os);
    return os;
}


// ************************************************************************* //
