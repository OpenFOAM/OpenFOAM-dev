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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::TroeFallOffFunction::TroeFallOffFunction
(
    const scalar alpha,
    const scalar Tsss,
    const scalar Ts,
    const scalar Tss
)
:
    alpha_(alpha),
    Tsss_(Tsss),
    Ts_(Ts),
    Tss_(Tss)
{}


inline Foam::TroeFallOffFunction::TroeFallOffFunction(Istream& is)
:
    alpha_(readScalar(is.readBegin("TroeFallOffFunction(Istream&)"))),
    Tsss_(readScalar(is)),
    Ts_(readScalar(is)),
    Tss_(readScalar(is))
{
    is.readEnd("TroeFallOffFunction(Istream&)");
}


inline Foam::TroeFallOffFunction::TroeFallOffFunction(const dictionary& dict)
:
    alpha_(readScalar(dict.lookup("alpha"))),
    Tsss_(readScalar(dict.lookup("Tsss"))),
    Ts_(readScalar(dict.lookup("Ts"))),
    Tss_(readScalar(dict.lookup("Tss")))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::TroeFallOffFunction::operator()
(
    const scalar T,
    const scalar Pr
) const
{
    scalar logFcent = log10
    (
        max
        (
            (1 - alpha_)*exp(-T/Tsss_) + alpha_*exp(-T/Ts_) + exp(-Tss_/T),
            SMALL
        )
    );

    scalar c = -0.4 - 0.67*logFcent;
    static const scalar d = 0.14;
    scalar n = 0.75 - 1.27*logFcent;

    scalar logPr = log10(max(Pr, SMALL));
    return pow(10.0, logFcent/(1.0 + sqr((logPr + c)/(n - d*(logPr + c)))));
}


inline void Foam::TroeFallOffFunction::write(Ostream& os) const
{
    os.writeKeyword("alpha") << alpha_ << token::END_STATEMENT << nl;
    os.writeKeyword("Tsss") << Tsss_ << token::END_STATEMENT << nl;
    os.writeKeyword("Ts") << Ts_ << token::END_STATEMENT << nl;
    os.writeKeyword("Tss") << Tss_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

inline Foam::Ostream& Foam::operator<<
(
    Foam::Ostream& os,
    const Foam::TroeFallOffFunction& tfof
)
{
    os  << token::BEGIN_LIST
        << tfof.alpha_
        << token::SPACE << tfof.Tsss_
        << token::SPACE << tfof.Ts_
        << token::SPACE << tfof.Tss_
        << token::END_LIST;

    return os;
}


// ************************************************************************* //
