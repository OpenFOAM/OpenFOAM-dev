/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2022 OpenFOAM Foundation
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

#include "MPLICcellStorage.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::pointField& Foam::MPLICcellStorage::points() const
{
    return points_;
}


inline const Foam::faceList& Foam::MPLICcellStorage::faces() const
{
    return faces_;
}


inline const Foam::edgeList& Foam::MPLICcellStorage::edges() const
{
    return edges_;
}


inline const Foam::labelListList& Foam::MPLICcellStorage::faceEdges() const
{
    return faceEdges_;
}

inline const Foam::labelList& Foam::MPLICcellStorage::cellPoints() const
{
    return cPoints_;
}


inline const Foam::labelList& Foam::MPLICcellStorage::cellFaces() const
{
    return cFaces_;
}


inline const Foam::labelList& Foam::MPLICcellStorage::cellEdges() const
{
    return cEdges_;
}


inline const Foam::boolList& Foam::MPLICcellStorage::isOwner() const
{
    return owns_;
}


inline const Foam::scalarField& Foam::MPLICcellStorage::pointsAlpha() const
{
    return pointsAlpha_;
}


inline const Foam::vectorField& Foam::MPLICcellStorage::pointsU() const
{
    return pointsU_;
}


inline Foam::scalar Foam::MPLICcellStorage::cellAlpha() const
{
    return cellAlpha_;
}


inline const Foam::vector& Foam::MPLICcellStorage::cellU() const
{
    return celllU_;
}


inline Foam::scalar Foam::MPLICcellStorage::V() const
{
    return volume_;
}


inline const Foam::vector& Foam::MPLICcellStorage::C() const
{
    return centre_;
}


inline const Foam::UIndirectList<Foam::vector>&
Foam::MPLICcellStorage::Sf() const
{
    return Sf_;
}


inline const Foam::UIndirectList<Foam::vector>&
Foam::MPLICcellStorage::Cf() const
{
    return Cf_;
}


inline const Foam::UIndirectList<Foam::scalar>&
Foam::MPLICcellStorage::magSf() const
{
    return magSf_;
}


inline Foam::scalar Foam::MPLICcellStorage::cellAlphaMax() const
{
    return cellAlphaMax_;
}


inline Foam::scalar Foam::MPLICcellStorage::cellAlphaMin() const
{
    return cellAlphaMin_;
}


inline const Foam::scalarField& Foam::MPLICcellStorage::facesAlphaMin() const
{
    return facesAlphaMin_;
}


inline const Foam::scalarField& Foam::MPLICcellStorage::facesAlphaMax() const
{
    return facesAlphaMax_;
}


inline Foam::label Foam::MPLICcellStorage::size() const
{
    return cFaces_.size();
}


// ************************************************************************* //
