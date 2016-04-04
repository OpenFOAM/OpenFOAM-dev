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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::label Foam::RBD::rigidBodyModel::nBodies() const
{
    return bodies_.size();
}


inline Foam::PtrList<Foam::RBD::rigidBody>
Foam::RBD::rigidBodyModel::bodies() const
{
    return bodies_;
}


inline const Foam::PtrList<Foam::RBD::joint>&
Foam::RBD::rigidBodyModel::joints() const
{
    return joints_;
}


inline Foam::label Foam::RBD::rigidBodyModel::nDoF() const
{
    return nDoF_;
}


inline Foam::label Foam::RBD::rigidBodyModel::nw() const
{
    return nw_;
}


inline const Foam::vector& Foam::RBD::rigidBodyModel::g() const
{
    return g_;
}


inline Foam::vector& Foam::RBD::rigidBodyModel::g()
{
    return g_;
}


inline const Foam::DynamicList<Foam::label>&
Foam::RBD::rigidBodyModel::lambda() const
{
    return lambda_;
}


inline const Foam::RBD::rigidBodyInertia&
Foam::RBD::rigidBodyModel::I(const label i) const
{
    return bodies_[i];
}


inline bool Foam::RBD::rigidBodyModel::merged(label bodyID) const
{
    return bodyID < 0;
}


inline Foam::label
Foam::RBD::rigidBodyModel::mergedBodyID(const label mergedBodyIndex) const
{
    return -1 - mergedBodyIndex;
}


inline Foam::label
Foam::RBD::rigidBodyModel::mergedBodyIndex(const label mergedBodyID) const
{
    return -1 - mergedBodyID;
}


inline const Foam::RBD::subBody&
Foam::RBD::rigidBodyModel::mergedBody(label mergedBodyID) const
{
    if (!merged(mergedBodyID))
    {
        FatalErrorInFunction
            << "Body " << mergedBodyID << " has not been merged"
            << abort(FatalError);
    }

    return mergedBodies_[mergedBodyIndex(mergedBodyID)];
}


inline Foam::label Foam::RBD::rigidBodyModel::bodyID(const word& name) const
{
    return bodyIDs_[name];
}


// ************************************************************************* //
