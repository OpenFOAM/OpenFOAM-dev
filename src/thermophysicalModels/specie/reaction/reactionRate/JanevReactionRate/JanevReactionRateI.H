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

inline Foam::JanevReactionRate::JanevReactionRate
(
    const scalar A,
    const scalar beta,
    const scalar Ta,
    const FixedList<scalar, nb_> b
)
:
    A_(A),
    beta_(beta),
    Ta_(Ta),
    b_(b)
{}


inline Foam::JanevReactionRate::JanevReactionRate
(
    const speciesTable&,
    Istream& is
)
:
    A_(readScalar(is.readBegin("JanevReactionRate(Istream&)"))),
    beta_(readScalar(is)),
    Ta_(readScalar(is)),
    b_(is)
{
    is.readEnd("JanevReactionRate(Istream&)");
}


inline Foam::JanevReactionRate::JanevReactionRate
(
    const speciesTable&,
    const dictionary& dict
)
:
    A_(readScalar(dict.lookup("A"))),
    beta_(readScalar(dict.lookup("beta"))),
    Ta_(readScalar(dict.lookup("Ta"))),
    b_(dict.lookup("b"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::JanevReactionRate::operator()
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

    scalar lnT = log(T);

    for (int n=0; n<nb_; n++)
    {
        expArg += b_[n]*pow(lnT, n);
    }

    lta *= exp(expArg);

    return lta;
}


inline void Foam::JanevReactionRate::write(Ostream& os) const
{
    os.writeKeyword("A") << A_ << nl;
    os.writeKeyword("beta") << beta_ << nl;
    os.writeKeyword("Ta") << Ta_ << nl;
    os.writeKeyword("b") << b_ << nl;
}


inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const JanevReactionRate& jrr
)
{
    os  << token::BEGIN_LIST
        << jrr.A_ << token::SPACE << jrr.beta_ << token::SPACE << jrr.Ta_;

    for (int n=0; n<JanevReactionRate::nb_; n++)
    {
        os  << token::SPACE << jrr.b_[n];
    }

    os << token::END_LIST;

    return os;
}


// ************************************************************************* //
