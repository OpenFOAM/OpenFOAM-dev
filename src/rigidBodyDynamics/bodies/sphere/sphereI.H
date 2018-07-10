/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline Foam::symmTensor Foam::RBD::sphere::I
(
    const scalar m,
    const scalar r
) const
{
    return ((2.0/5.0)*m*sqr(r))*Foam::I;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::RBD::sphere::sphere
(
    const word& name,
    const scalar m,
    const vector& c,
    const scalar r
)
:
    rigidBody(name, m, c, I(m, r)),
    r_(r)
{}


inline Foam::RBD::sphere::sphere
(
    const word& name,
    const dictionary& dict
)
:
    rigidBody(name, rigidBodyInertia()),
    r_(readScalar(dict.lookup("radius")))
{
    const scalar m(readScalar(dict.lookup("mass")));
    const vector c(dict.lookup("centreOfMass"));
    rigidBodyInertia::operator=(rigidBodyInertia(m, c, I(m, r_)));
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline Foam::scalar Foam::RBD::sphere::r() const
{
    return r_;
}


// ************************************************************************* //
