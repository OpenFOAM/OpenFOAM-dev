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

inline Foam::LindemannFallOffFunction::LindemannFallOffFunction()
{}


inline Foam::LindemannFallOffFunction::LindemannFallOffFunction
(
    const dictionary& dict
)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::LindemannFallOffFunction::operator()
(
    const scalar,
    const scalar
) const
{
    return 1.0;
}


inline Foam::scalar Foam::LindemannFallOffFunction::ddT
(
    const scalar,
    const scalar,
    const scalar,
    const scalar
) const
{
    return 0;
}


inline Foam::scalar Foam::LindemannFallOffFunction::ddc
(
    const scalar,
    const scalar,
    const scalar,
    const scalar
) const
{
    return 0;
}


inline void Foam::LindemannFallOffFunction::write(Ostream& os) const
{}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

inline Foam::Ostream& Foam::operator<<
(
    Foam::Ostream& os,
    const Foam::LindemannFallOffFunction&
)
{
    return os;
}


// ************************************************************************* //
