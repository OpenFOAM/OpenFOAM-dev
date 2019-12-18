/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

inline Foam::thirdBodyArrheniusReactionRate::thirdBodyArrheniusReactionRate
(
    const scalar A,
    const scalar beta,
    const scalar Ta,
    const thirdBodyEfficiencies& tbes
)
:
    ArrheniusReactionRate(A, beta, Ta),
    thirdBodyEfficiencies_(tbes)
{
    forAll(tbes, i)
    {
        beta_.append(Tuple2<label, scalar>(i, tbes[i]));
    }
}


inline Foam::thirdBodyArrheniusReactionRate::thirdBodyArrheniusReactionRate
(
    const speciesTable& species,
    const dictionary& dict
)
:
    ArrheniusReactionRate
    (
        species,
        dict
    ),
    thirdBodyEfficiencies_(species, dict)
{
    forAll(thirdBodyEfficiencies_, i)
    {
        beta_.append
        (
            Tuple2<label, scalar>
            (
                i,
                thirdBodyEfficiencies_[i]
            )
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::thirdBodyArrheniusReactionRate::operator()
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return
        thirdBodyEfficiencies_.M(c)
       *ArrheniusReactionRate::operator()(p, T, c, li);
}


inline Foam::scalar Foam::thirdBodyArrheniusReactionRate::ddT
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return
        thirdBodyEfficiencies_.M(c)
       *ArrheniusReactionRate::ddT(p, T, c, li);
}


inline void Foam::thirdBodyArrheniusReactionRate::dcidc
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    scalarField& dcidc
) const
{
    scalar M = thirdBodyEfficiencies_.M(c);
    forAll(beta_, i)
    {
        dcidc[i] = beta_[i].second()/max(M, small);
    }
}


inline Foam::scalar Foam::thirdBodyArrheniusReactionRate::dcidT
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return -1.0/T;
}


inline void Foam::thirdBodyArrheniusReactionRate::write(Ostream& os) const
{
    ArrheniusReactionRate::write(os);
    thirdBodyEfficiencies_.write(os);
}


inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const thirdBodyArrheniusReactionRate& arr
)
{
    arr.write(os);
    return os;
}


// ************************************************************************* //
