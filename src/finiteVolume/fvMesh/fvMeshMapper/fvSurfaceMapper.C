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

Description
    FV surface mapper.

\*---------------------------------------------------------------------------*/

#include "fvSurfaceMapper.H"
#include "fvMesh.H"
#include "mapPolyMesh.H"
#include "faceMapper.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fvSurfaceMapper::calcAddressing() const
{
    if
    (
        directAddrPtr_
     || interpolationAddrPtr_
     || weightsPtr_
     || insertedObjectLabelsPtr_
    )
    {
        FatalErrorIn("void fvSurfaceMapper::calcAddressing() const)")
            << "Addressing already calculated"
            << abort(FatalError);
    }

    // Mapping

    const label oldNInternal = faceMap_.nOldInternalFaces();

    // Assemble the maps
    if (direct())
    {
        // Direct mapping - slice to size
        directAddrPtr_ =
            new labelList
            (
                labelList::subList(faceMap_.directAddressing(), size())
            );
        labelList& addr = *directAddrPtr_;

        // Adjust for creation of an internal face from a boundary face
        forAll(addr, faceI)
        {
            if (addr[faceI] > oldNInternal)
            {
                addr[faceI] = 0;
            }
        }
    }
    else
    {
        // Interpolative mapping - slice to size
        interpolationAddrPtr_ =
            new labelListList
            (
                labelListList::subList(faceMap_.addressing(), size())
            );
        labelListList& addr = *interpolationAddrPtr_;

        weightsPtr_ =
            new scalarListList
            (
                scalarListList::subList(faceMap_.weights(), size())
            );
        scalarListList& w = *weightsPtr_;

        // Adjust for creation of an internal face from a boundary face
        forAll(addr, faceI)
        {
            if (max(addr[faceI]) >= oldNInternal)
            {
                addr[faceI] = labelList(1, label(0));
                w[faceI] = scalarList(1, 1.0);
            }
        }
    }

    // Inserted objects

    // If there are, assemble the labels
    if (insertedObjects())
    {
        const labelList& insFaces = faceMap_.insertedObjectLabels();

        insertedObjectLabelsPtr_ = new labelList(insFaces.size());
        labelList& ins = *insertedObjectLabelsPtr_;

        label nIns = 0;

        forAll(insFaces, faceI)
        {
            // If the face is internal, keep it here
            if (insFaces[faceI] < size())
            {
                ins[nIns] = insFaces[faceI];
                nIns++;
            }
        }

        ins.setSize(nIns);
    }
    else
    {
        // No inserted objects
        insertedObjectLabelsPtr_ = new labelList(0);
    }
}


void Foam::fvSurfaceMapper::clearOut()
{
    deleteDemandDrivenData(directAddrPtr_);
    deleteDemandDrivenData(interpolationAddrPtr_);
    deleteDemandDrivenData(weightsPtr_);

    deleteDemandDrivenData(insertedObjectLabelsPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::fvSurfaceMapper::fvSurfaceMapper
(
    const fvMesh& mesh,
    const faceMapper& fMapper
)
:
    mesh_(mesh),
    faceMap_(fMapper),
    directAddrPtr_(NULL),
    interpolationAddrPtr_(NULL),
    weightsPtr_(NULL),
    insertedObjectLabelsPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvSurfaceMapper::~fvSurfaceMapper()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelUList& Foam::fvSurfaceMapper::directAddressing() const
{
    if (!direct())
    {
        FatalErrorIn
        (
            "const labelUList& fvSurfaceMapper::"
            "directAddressing() const"
        )   << "Requested direct addressing for an interpolative mapper."
            << abort(FatalError);
    }

    if (!directAddrPtr_)
    {
        calcAddressing();
    }

    return *directAddrPtr_;
}


const Foam::labelListList& Foam::fvSurfaceMapper::addressing() const
{
    if (direct())
    {
        FatalErrorIn
        (
            "const labelListList& fvSurfaceMapper::addressing() const"
        )   << "Requested interpolative addressing for a direct mapper."
            << abort(FatalError);
    }

    if (!interpolationAddrPtr_)
    {
        calcAddressing();
    }

    return *interpolationAddrPtr_;
}


const Foam::scalarListList& Foam::fvSurfaceMapper::weights() const
{
    if (direct())
    {
        FatalErrorIn
        (
            "const scalarListList& fvSurfaceMapper::weights() const"
        )   << "Requested interpolative weights for a direct mapper."
            << abort(FatalError);
    }

    if (!weightsPtr_)
    {
        calcAddressing();
    }

    return *weightsPtr_;
}


const Foam::labelList& Foam::fvSurfaceMapper::insertedObjectLabels() const
{
    if (!insertedObjectLabelsPtr_)
    {
        calcAddressing();
    }

    return *insertedObjectLabelsPtr_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// ************************************************************************* //
