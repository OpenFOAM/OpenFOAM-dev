/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
#include "polyTopoChangeMap.H"
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
        FatalErrorInFunction
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
                labelList::subList
                (
                    faceMap_.directAddressing(),
                    mesh_.nInternalFaces()
                )
            );
        labelList& addr = *directAddrPtr_;

        // Adjust for creation of an internal face from a boundary face
        forAll(addr, facei)
        {
            if (addr[facei] > oldNInternal)
            {
                addr[facei] = 0;
            }
        }
    }
    else
    {
        // Interpolative mapping - slice to size
        interpolationAddrPtr_ =
            new labelListList
            (
                labelListList::subList
                (
                    faceMap_.addressing(),
                    mesh_.nInternalFaces()
                )
            );
        labelListList& addr = *interpolationAddrPtr_;

        weightsPtr_ =
            new scalarListList
            (
                scalarListList::subList
                (
                    faceMap_.weights(),
                    mesh_.nInternalFaces()
                )
            );
        scalarListList& w = *weightsPtr_;

        // Adjust for creation of an internal face from a boundary face
        forAll(addr, facei)
        {
            if (max(addr[facei]) >= oldNInternal)
            {
                addr[facei] = labelList(1, label(0));
                w[facei] = scalarList(1, 1.0);
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

        forAll(insFaces, facei)
        {
            // If the face is internal, keep it here
            if (insFaces[facei] < mesh_.nInternalFaces())
            {
                ins[nIns] = insFaces[facei];
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

Foam::fvSurfaceMapper::fvSurfaceMapper
(
    const fvMesh& mesh,
    const faceMapper& fMapper
)
:
    mesh_(mesh),
    faceMap_(fMapper),
    directAddrPtr_(nullptr),
    interpolationAddrPtr_(nullptr),
    weightsPtr_(nullptr),
    insertedObjectLabelsPtr_(nullptr)
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
        FatalErrorInFunction
            << "Requested direct addressing for an interpolative mapper."
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
        FatalErrorInFunction
            << "Requested interpolative addressing for a direct mapper."
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
        FatalErrorInFunction
            << "Requested interpolative weights for a direct mapper."
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
