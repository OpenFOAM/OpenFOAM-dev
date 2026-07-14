/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2026 OpenFOAM Foundation
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

#include "logPolynomialTransport.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo, int PolySize>
Foam::logPolynomialTransport<Thermo, PolySize>::logPolynomialTransport
(
    const word& name,
    const dictionary& dict,
    const dictionary& subDict
)
:
    Thermo(name, dict),
    muLogCoeffs_
    (
        subDict.lookup<FixedPolynomial<scalar, PolySize>>
        (
            "muLogCoeffs<" + Foam::name(PolySize) + '>',
            Function1s::unitSets({units::none, units::none})
        )
    ),
    kappaLogCoeffs_
    (
        subDict.lookup<FixedPolynomial<scalar, PolySize>>
        (
            "kappaLogCoeffs<" + Foam::name(PolySize) + '>',
            Function1s::unitSets({units::none, units::none})
        )
    )
{}


template<class Thermo, int PolySize>
Foam::logPolynomialTransport<Thermo, PolySize>::logPolynomialTransport
(
    const word& name,
    const dictionary& dict
)
:
    logPolynomialTransport(name, dict, dict.subDict("transport"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo, int PolySize>
void Foam::logPolynomialTransport<Thermo, PolySize>::write(Ostream& os) const
{
    Thermo::write(os);

    writeEntry
    (
        os,
        "transport",
        dictionary::entries
        (
            word("muLogCoeffs<" + Foam::name(PolySize) + '>'),
            muLogCoeffs_,
            word("kappaLogCoeffs<" + Foam::name(PolySize) + '>'),
            kappaLogCoeffs_
        )
    );
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Thermo, int PolySize>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const logPolynomialTransport<Thermo, PolySize>& pt
)
{
    pt.write(os);
    return os;
}


// ************************************************************************* //
