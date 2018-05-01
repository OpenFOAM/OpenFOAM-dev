/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

inline Foam::solidArrheniusReactionRate::solidArrheniusReactionRate
(
    const scalar A,
    const scalar Ta,
    const scalar Tcrit
    // const scalar nReact
)
:
    A_(A),
    Ta_(Ta),
    Tcrit_(Tcrit)
{}


inline Foam::solidArrheniusReactionRate::solidArrheniusReactionRate
(
    const speciesTable&,
    const dictionary& dict
)
:
    A_(readScalar(dict.lookup("A"))),
    Ta_(readScalar(dict.lookup("Ta"))),
    Tcrit_(readScalar(dict.lookup("Tcrit")))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::solidArrheniusReactionRate::operator()
(
    const scalar,
    const scalar T,
    const scalarField&
) const
{
    scalar ak = A_;

    if (T < Tcrit_)
    {
        ak *= 0.0;
    }
    else
    {
        ak *= exp(-Ta_/T);
    }

    return ak;
}


inline void Foam::solidArrheniusReactionRate::write(Ostream& os) const
{
    os.writeKeyword("A") << A_ << token::END_STATEMENT << nl;
    os.writeKeyword("Ta") << Ta_ << token::END_STATEMENT << nl;
    os.writeKeyword("Tcrit") << Tcrit_ << token::END_STATEMENT << nl;
}


inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const solidArrheniusReactionRate& arr
)
{
    arr.write(os);
    return os;
}


// ************************************************************************* //
