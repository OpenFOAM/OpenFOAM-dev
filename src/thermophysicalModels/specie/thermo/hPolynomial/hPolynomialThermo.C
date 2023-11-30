/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "hPolynomialThermo.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState, int PolySize>
Foam::hPolynomialThermo<EquationOfState, PolySize>::hPolynomialThermo
(
    const word& name,
    const dictionary& dict
)
:
    EquationOfState(name, dict),
    hf_
    (
        dict
       .subDict("thermodynamics")
       .lookupBackwardsCompatible<scalar>({"hf", "Hf"})
    ),
    sf_
    (
        dict
       .subDict("thermodynamics")
       .lookupBackwardsCompatible<scalar>({"sf", "Sf"})
    ),
    CpCoeffs_
    (
        dict.subDict("thermodynamics").lookup
        (
            "CpCoeffs<" + Foam::name(PolySize) + '>'
        )
    ),
    hCoeffs_(),
    sCoeffs_()
{
    hCoeffs_ = CpCoeffs_.integral();
    sCoeffs_ = CpCoeffs_.integralMinus1();

    // Offset h poly so that it is relative to the enthalpy at Tstd
    hCoeffs_[0] += hf_ - hCoeffs_.value(Tstd);

    // Offset s poly so that it is relative to the entropy at Tstd
    sCoeffs_[0] += sf_ - sCoeffs_.value(Tstd);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState, int PolySize>
void Foam::hPolynomialThermo<EquationOfState, PolySize>::write
(
    Ostream& os
) const
{
    EquationOfState::write(os);

    dictionary dict("thermodynamics");
    dict.add("hf", hf_);
    dict.add("sf", sf_);
    dict.add
    (
        word("CpCoeffs<" + Foam::name(PolySize) + '>'),
        CpCoeffs_
    );
    os  << indent << dict.dictName() << dict;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class EquationOfState, int PolySize>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const hPolynomialThermo<EquationOfState, PolySize>& pt
)
{
    pt.write(os);
    return os;
}


// ************************************************************************* //
