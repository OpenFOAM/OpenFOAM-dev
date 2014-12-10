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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::LandauTellerReactionRate::LandauTellerReactionRate
(
    const scalar A,
    const scalar beta,
    const scalar Ta,
    const scalar B,
    const scalar C
)
:
    A_(A),
    beta_(beta),
    Ta_(Ta),
    B_(B),
    C_(C)
{}


inline Foam::LandauTellerReactionRate::LandauTellerReactionRate
(
    const speciesTable&,
    Istream& is
)
:
    A_(readScalar(is.readBegin("LandauTellerReactionRate(Istream&)"))),
    beta_(readScalar(is)),
    Ta_(readScalar(is)),
    B_(readScalar(is)),
    C_(readScalar(is))
{
    is.readEnd("LandauTellerReactionRate(Istream&)");
}


inline Foam::LandauTellerReactionRate::LandauTellerReactionRate
(
    const speciesTable&,
    const dictionary& dict
)
:
    A_(readScalar(dict.lookup("A"))),
    beta_(readScalar(dict.lookup("beta"))),
    Ta_(readScalar(dict.lookup("Ta"))),
    B_(readScalar(dict.lookup("B"))),
    C_(readScalar(dict.lookup("C")))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::LandauTellerReactionRate::operator()
(
    const scalar p,
    const scalar T,
    const scalarField&
) const
{
    scalar lta = A_;

    if (mag(beta_) > VSMALL)
    {
        lta *= pow(T, beta_);
    }

    scalar expArg = 0.0;

    if (mag(Ta_) > VSMALL)
    {
        expArg -= Ta_/T;
    }

    if (mag(B_) > VSMALL)
    {
        expArg += B_/cbrt(T);
    }

    if (mag(C_) > VSMALL)
    {
        expArg += C_/pow(T, 2.0/3.0);
    }

    if (mag(expArg) > VSMALL)
    {
        lta *= exp(expArg);
    }

    return lta;
}


inline void Foam::LandauTellerReactionRate::write(Ostream& os) const
{
    os.writeKeyword("A") << A_ << token::END_STATEMENT << nl;
    os.writeKeyword("beta") << beta_ << token::END_STATEMENT << nl;
    os.writeKeyword("Ta") << Ta_ << token::END_STATEMENT << nl;
    os.writeKeyword("B") << B_ << token::END_STATEMENT << nl;
    os.writeKeyword("C") << C_ << token::END_STATEMENT << nl;
}


inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const LandauTellerReactionRate& arr
)
{
    os  << token::BEGIN_LIST
        << arr.A_ << token::SPACE << arr.beta_ << token::SPACE << arr.Ta_
        << token::SPACE << arr.B_ << token::SPACE << arr.C_
        << token::END_LIST;
    return os;
}


// ************************************************************************* //
