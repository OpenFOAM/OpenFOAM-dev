/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2024 OpenFOAM Foundation
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
    const autoPtr<rigidBody>& bodyPtr,
    const word& masterName,
    const label masterID,
    const spatialTransform& masterXT
)
:
    body_(bodyPtr),
    masterName_(masterName),
    masterIndex_(masterID),
    masterXT_(masterXT)
{}


inline Foam::autoPtr<Foam::RBD::subBody> Foam::RBD::subBody::clone() const
{
    return autoPtr<subBody>(new subBody(*this));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::RBD::rigidBody& Foam::RBD::subBody::body() const
{
    return body_();
}


inline const Foam::word& Foam::RBD::subBody::name() const
{
    return body_->name();
}


inline const Foam::word& Foam::RBD::subBody::masterName() const
{
    return masterName_;
}


inline Foam::label Foam::RBD::subBody::masterIndex() const
{
    return masterIndex_;
}


inline const Foam::spatialTransform& Foam::RBD::subBody::masterXT() const
{
    return masterXT_;
}


// ************************************************************************* //
