/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Istream& is
)
:
    EquationOfState(is),
    Hf_(readScalar(is)),
    Sf_(readScalar(is)),
    CpCoeffs_("CpCoeffs<" + Foam::name(PolySize) + '>', is),
    hCoeffs_(),
    sCoeffs_()
{
    Hf_ *= this->W();
    Sf_ *= this->W();
    CpCoeffs_ *= this->W();

    hCoeffs_ = CpCoeffs_.integral();
    sCoeffs_ = CpCoeffs_.integralMinus1();

    // Offset h poly so that it is relative to the enthalpy at Tstd
    hCoeffs_[0] += Hf_ - hCoeffs_.value(specie::Tstd);

    // Offset s poly so that it is relative to the entropy at Tstd
    sCoeffs_[0] += Sf_ - sCoeffs_.value(specie::Tstd);
}


template<class EquationOfState, int PolySize>
Foam::hPolynomialThermo<EquationOfState, PolySize>::hPolynomialThermo
(
    const dictionary& dict
)
:
    EquationOfState(dict),
    Hf_(readScalar(dict.subDict("thermodynamics").lookup("Hf"))),
    Sf_(readScalar(dict.subDict("thermodynamics").lookup("Sf"))),
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
    Hf_ *= this->W();
    Sf_ *= this->W();
    CpCoeffs_ *= this->W();

    hCoeffs_ = CpCoeffs_.integral();
    sCoeffs_ = CpCoeffs_.integralMinus1();

    // Offset h poly so that it is relative to the enthalpy at Tstd
    hCoeffs_[0] += Hf_ - hCoeffs_.value(specie::Tstd);

    // Offset s poly so that it is relative to the entropy at Tstd
    sCoeffs_[0] += Sf_ - sCoeffs_.value(specie::Tstd);
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
    dict.add("Hf", Hf_/this->W());
    dict.add("Sf", Sf_/this->W());
    dict.add
    (
        word("CpCoeffs<" + Foam::name(PolySize) + '>'),
        CpCoeffs_/this->W()
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
    os  << static_cast<const EquationOfState&>(pt) << tab
        << pt.Hf_/pt.W() << tab
        << pt.Sf_/pt.W() << tab
        << "CpCoeffs<" << Foam::name(PolySize) << '>' << tab
        << pt.CpCoeffs_/pt.W();

    os.check
    (
        "operator<<"
        "("
            "Ostream&, "
            "const hPolynomialThermo<EquationOfState, PolySize>&"
        ")"
    );

    return os;
}


// ************************************************************************* //
