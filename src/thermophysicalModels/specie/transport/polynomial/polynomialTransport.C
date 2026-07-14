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

#include "polynomialTransport.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo, int PolySize>
Foam::polynomialTransport<Thermo, PolySize>::polynomialTransport
(
    const word& name,
    const dictionary& dict,
    const dictionary& subDict
)
:
    Thermo(name, dict),
    muCoeffs_
    (
        subDict.lookup<FixedPolynomial<scalar, PolySize>>
        (
            "muCoeffs<" + Foam::name(PolySize) + '>',
            Function1s::unitSets({dimTemperature, dimDynamicViscosity})
        )
    ),
    kappaCoeffs_
    (
        subDict.lookup<FixedPolynomial<scalar, PolySize>>
        (
            "kappaCoeffs<" + Foam::name(PolySize) + '>',
            Function1s::unitSets({dimTemperature, dimThermalConductivity})
        )
    )
{}


template<class Thermo, int PolySize>
Foam::polynomialTransport<Thermo, PolySize>::polynomialTransport
(
    const word& name,
    const dictionary& dict
)
:
    polynomialTransport(name, dict, dict.subDict("transport"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo, int PolySize>
void Foam::polynomialTransport<Thermo, PolySize>::write(Ostream& os) const
{
    Thermo::write(os);

    writeEntry
    (
        os,
        "transport",
        dictionary::entries
        (
            word("muCoeffs<" + Foam::name(PolySize) + '>'), muCoeffs_,
            word("kappaCoeffs<" + Foam::name(PolySize) + '>'), kappaCoeffs_
        )
    );
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Thermo, int PolySize>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const polynomialTransport<Thermo, PolySize>& pt
)
{
    pt.write(os);
    return os;
}


// ************************************************************************* //
