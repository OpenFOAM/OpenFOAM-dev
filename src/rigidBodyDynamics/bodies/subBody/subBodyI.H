/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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

inline Foam::RBD::subBody::subBody
(
    autoPtr<rigidBody> bodyPtr,
    const label parentID,
    const spatialTransform& parentXT
)
:
    rigidBody(bodyPtr()),
    body_(bodyPtr),
    parentID_(parentID),
    parentXT_(parentXT)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::RBD::rigidBody& Foam::RBD::subBody::body() const
{
    return body_();
}


inline Foam::label Foam::RBD::subBody::parentID() const
{
    return parentID_;
}


inline const Foam::spatialTransform& Foam::RBD::subBody::parentXT() const
{
    return parentXT_;
}


// ************************************************************************* //
