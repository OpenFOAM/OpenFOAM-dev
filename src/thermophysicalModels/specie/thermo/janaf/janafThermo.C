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

#include "janafThermo.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class EquationOfState>
void Foam::janafThermo<EquationOfState>::checkInputData() const
{
    if (Tlow_ >= Thigh_)
    {
        FatalErrorIn("janafThermo<EquationOfState>::check()")
            << "Tlow(" << Tlow_ << ") >= Thigh(" << Thigh_ << ')'
            << exit(FatalError);
    }

    if (Tcommon_ <= Tlow_)
    {
        FatalErrorIn("janafThermo<EquationOfState>::check()")
            << "Tcommon(" << Tcommon_ << ") <= Tlow(" << Tlow_ << ')'
            << exit(FatalError);
    }

    if (Tcommon_ > Thigh_)
    {
        FatalErrorIn("janafThermo<EquationOfState>::check()")
            << "Tcommon(" << Tcommon_ << ") > Thigh(" << Thigh_ << ')'
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::janafThermo<EquationOfState>::janafThermo(Istream& is)
:
    EquationOfState(is),
    Tlow_(readScalar(is)),
    Thigh_(readScalar(is)),
    Tcommon_(readScalar(is))
{
    checkInputData();

    forAll(highCpCoeffs_, i)
    {
        is >> highCpCoeffs_[i];
    }

    forAll(lowCpCoeffs_, i)
    {
        is >> lowCpCoeffs_[i];
    }

    // Check state of Istream
    is.check("janafThermo::janafThermo(Istream& is)");
}


template<class EquationOfState>
Foam::janafThermo<EquationOfState>::janafThermo(const dictionary& dict)
:
    EquationOfState(dict),
    Tlow_(readScalar(dict.subDict("thermodynamics").lookup("Tlow"))),
    Thigh_(readScalar(dict.subDict("thermodynamics").lookup("Thigh"))),
    Tcommon_(readScalar(dict.subDict("thermodynamics").lookup("Tcommon"))),
    highCpCoeffs_(dict.subDict("thermodynamics").lookup("highCpCoeffs")),
    lowCpCoeffs_(dict.subDict("thermodynamics").lookup("lowCpCoeffs"))
{
    checkInputData();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
void Foam::janafThermo<EquationOfState>::write(Ostream& os) const
{
    EquationOfState::write(os);

    dictionary dict("thermodynamics");
    dict.add("Tlow", Tlow_);
    dict.add("Thigh", Thigh_);
    dict.add("Tcommon", Tcommon_);
    dict.add("highCpCoeffs", highCpCoeffs_);
    dict.add("lowCpCoeffs", lowCpCoeffs_);
    os  << indent << dict.dictName() << dict;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const janafThermo<EquationOfState>& jt
)
{
    os  << static_cast<const EquationOfState&>(jt) << nl
        << "    " << jt.Tlow_
        << tab << jt.Thigh_
        << tab << jt.Tcommon_;

    os << nl << "    ";

    forAll(jt.highCpCoeffs_, i)
    {
        os << jt.highCpCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(jt.lowCpCoeffs_, i)
    {
        os << jt.lowCpCoeffs_[i] << ' ';
    }

    os << endl;

    os.check
    (
        "operator<<(Ostream& os, const janafThermo<EquationOfState>& jt)"
    );

    return os;
}


// ************************************************************************* //
