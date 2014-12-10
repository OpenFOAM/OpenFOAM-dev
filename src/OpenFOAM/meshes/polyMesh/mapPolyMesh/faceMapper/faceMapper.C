/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "faceMapper.H"
#include "demandDrivenData.H"
#include "polyMesh.H"
#include "mapPolyMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::faceMapper::calcAddressing() const
{
    if
    (
        directAddrPtr_
     || interpolationAddrPtr_
     || weightsPtr_
     || insertedFaceLabelsPtr_
    )
    {
        FatalErrorIn("void faceMapper::calcAddressing() const")
            << "Addressing already calculated."
            << abort(FatalError);
    }

    if (direct())
    {
        // Direct addressing, no weights

        directAddrPtr_ = new labelList(mpm_.faceMap());
        labelList& directAddr = *directAddrPtr_;

        // Reset the size of addressing list to contain only live faces
        directAddr.setSize(mesh_.nFaces());

        insertedFaceLabelsPtr_ = new labelList(mesh_.nFaces());
        labelList& insertedFaces = *insertedFaceLabelsPtr_;

        label nInsertedFaces = 0;

        forAll(directAddr, faceI)
        {
            if (directAddr[faceI] < 0)
            {
                // Found inserted face
                directAddr[faceI] = 0;
                insertedFaces[nInsertedFaces] = faceI;
                nInsertedFaces++;
            }
        }

        insertedFaces.setSize(nInsertedFaces);
    }
    else
    {
        // Interpolative addressing

        interpolationAddrPtr_ = new labelListList(mesh_.nFaces());
        labelListList& addr = *interpolationAddrPtr_;

        weightsPtr_ = new scalarListList(mesh_.nFaces());
        scalarListList& w = *weightsPtr_;

        const List<objectMap>& ffp = mpm_.facesFromPointsMap();

        forAll(ffp, ffpI)
        {
            // Get addressing
            const labelList& mo = ffp[ffpI].masterObjects();

            label faceI = ffp[ffpI].index();

            if (addr[faceI].size())
            {
                FatalErrorIn("void faceMapper::calcAddressing() const")
                    << "Master face " << faceI
                    << " mapped from point faces " << mo
                    << " already destination of mapping." << abort(FatalError);
            }

            // Map from masters, uniform weights
            addr[faceI] = mo;
            w[faceI] = scalarList(mo.size(), 1.0/mo.size());
        }

        const List<objectMap>& ffe = mpm_.facesFromEdgesMap();

        forAll(ffe, ffeI)
        {
            // Get addressing
            const labelList& mo = ffe[ffeI].masterObjects();

            label faceI = ffe[ffeI].index();

            if (addr[faceI].size())
            {
                FatalErrorIn("void faceMapper::calcAddressing() const")
                    << "Master face " << faceI
                    << " mapped from edge faces " << mo
                    << " already destination of mapping." << abort(FatalError);
            }

            // Map from masters, uniform weights
            addr[faceI] = mo;
            w[faceI] = scalarList(mo.size(), 1.0/mo.size());
        }

        const List<objectMap>& fff = mpm_.facesFromFacesMap();

        forAll(fff, fffI)
        {
            // Get addressing
            const labelList& mo = fff[fffI].masterObjects();

            label faceI = fff[fffI].index();

            if (addr[faceI].size())
            {
                FatalErrorIn("void faceMapper::calcAddressing() const")
                    << "Master face " << faceI
                    << " mapped from face faces " << mo
                    << " already destination of mapping." << abort(FatalError);
            }

            // Map from masters, uniform weights
            addr[faceI] = mo;
            w[faceI] = scalarList(mo.size(), 1.0/mo.size());
        }


        // Do mapped faces. Note that can already be set from facesFromFaces
        // so check if addressing size still zero.
        const labelList& fm = mpm_.faceMap();

        forAll(fm, faceI)
        {
            if (fm[faceI] > -1 && addr[faceI].empty())
            {
                // Mapped from a single face
                addr[faceI] = labelList(1, fm[faceI]);
                w[faceI] = scalarList(1, 1.0);
            }
        }


        // Grab inserted faces (for them the size of addressing is still zero)

        insertedFaceLabelsPtr_ = new labelList(mesh_.nFaces());
        labelList& insertedFaces = *insertedFaceLabelsPtr_;

        label nInsertedFaces = 0;

        forAll(addr, faceI)
        {
            if (addr[faceI].empty())
            {
                // Mapped from a dummy face
                addr[faceI] = labelList(1, label(0));
                w[faceI] = scalarList(1, 1.0);

                insertedFaces[nInsertedFaces] = faceI;
                nInsertedFaces++;
            }
        }

        insertedFaces.setSize(nInsertedFaces);
    }
}


void Foam::faceMapper::clearOut()
{
    deleteDemandDrivenData(directAddrPtr_);
    deleteDemandDrivenData(interpolationAddrPtr_);
    deleteDemandDrivenData(weightsPtr_);
    deleteDemandDrivenData(insertedFaceLabelsPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::faceMapper::faceMapper(const mapPolyMesh& mpm)
:
    mesh_(mpm.mesh()),
    mpm_(mpm),
    insertedFaces_(true),
    direct_(false),
    directAddrPtr_(NULL),
    interpolationAddrPtr_(NULL),
    weightsPtr_(NULL),
    insertedFaceLabelsPtr_(NULL)
{
    // Check for possibility of direct mapping
    if
    (
        mpm_.facesFromPointsMap().empty()
     && mpm_.facesFromEdgesMap().empty()
     && mpm_.facesFromFacesMap().empty()
    )
    {
        direct_ = true;
    }
    else
    {
        direct_ = false;
    }

    // Check for inserted faces
    if (direct_ && (mpm_.faceMap().empty() || min(mpm_.faceMap()) > -1))
    {
        insertedFaces_ = false;
    }
    else
    {
        // Need to check all 3 lists to see if there are inserted faces
        // with no owner

        // Make a copy of the face map, add the entries for faces from points
        // and faces from edges and check for left-overs
        labelList fm(mesh_.nFaces(), -1);

        const List<objectMap>& ffp = mpm_.facesFromPointsMap();

        forAll(ffp, ffpI)
        {
            fm[ffp[ffpI].index()] = 0;
        }

        const List<objectMap>& ffe = mpm_.facesFromEdgesMap();

        forAll(ffe, ffeI)
        {
            fm[ffe[ffeI].index()] = 0;
        }

        const List<objectMap>& fff = mpm_.facesFromFacesMap();

        forAll(fff, fffI)
        {
            fm[fff[fffI].index()] = 0;
        }

        if (min(fm) < 0)
        {
            insertedFaces_ = true;
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::faceMapper::~faceMapper()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::faceMapper::size() const
{
    return mesh_.nFaces();
}


Foam::label Foam::faceMapper::sizeBeforeMapping() const
{
    return mpm_.nOldFaces();
}


Foam::label Foam::faceMapper::internalSizeBeforeMapping() const
{
    return mpm_.nOldInternalFaces();
}


const Foam::labelUList& Foam::faceMapper::directAddressing() const
{
    if (!direct())
    {
        FatalErrorIn
        (
            "const labelUList& faceMapper::directAddressing() const"
        )   << "Requested direct addressing for an interpolative mapper."
            << abort(FatalError);
    }

    if (!insertedObjects())
    {
        // No inserted faces.  Re-use faceMap
        return mpm_.faceMap();
    }
    else
    {
        if (!directAddrPtr_)
        {
            calcAddressing();
        }

        return *directAddrPtr_;
    }
}


const Foam::labelListList& Foam::faceMapper::addressing() const
{
    if (direct())
    {
        FatalErrorIn
        (
            "const labelListList& faceMapper::addressing() const"
        )   << "Requested interpolative addressing for a direct mapper."
            << abort(FatalError);
    }

    if (!interpolationAddrPtr_)
    {
        calcAddressing();
    }

    return *interpolationAddrPtr_;
}


const Foam::scalarListList& Foam::faceMapper::weights() const
{
    if (direct())
    {
        FatalErrorIn
        (
            "const scalarListList& faceMapper::weights() const"
        )   << "Requested interpolative weights for a direct mapper."
            << abort(FatalError);
    }

    if (!weightsPtr_)
    {
        calcAddressing();
    }

    return *weightsPtr_;
}


const Foam::labelList& Foam::faceMapper::insertedObjectLabels() const
{
    if (!insertedFaceLabelsPtr_)
    {
        if (!insertedObjects())
        {
            // There are no inserted faces
            insertedFaceLabelsPtr_ = new labelList(0);
        }
        else
        {
            calcAddressing();
        }
    }

    return *insertedFaceLabelsPtr_;
}


const Foam::labelHashSet& Foam::faceMapper::flipFaceFlux() const
{
    return mpm_.flipFaceFlux();
}


Foam::label Foam::faceMapper::nOldInternalFaces() const
{
    return mpm_.nOldInternalFaces();
}


const Foam::labelList& Foam::faceMapper::oldPatchStarts() const
{
    return mpm_.oldPatchStarts();
}


const Foam::labelList& Foam::faceMapper::oldPatchSizes() const
{
    return mpm_.oldPatchSizes();
}


// ************************************************************************* //
