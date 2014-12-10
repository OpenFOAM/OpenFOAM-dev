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

\*---------------------------------------------------------------------------*/

#include "pointMapper.H"
#include "demandDrivenData.H"
#include "pointMesh.H"
#include "mapPolyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::pointMapper::calcAddressing() const
{
    if
    (
        directAddrPtr_
     || interpolationAddrPtr_
     || weightsPtr_
     || insertedPointLabelsPtr_
    )
    {
        FatalErrorIn("void pointMapper::calcAddressing() const")
            << "Addressing already calculated."
            << abort(FatalError);
    }

    if (direct())
    {
        // Direct addressing, no weights

        directAddrPtr_ = new labelList(mpm_.pointMap());
        labelList& directAddr = *directAddrPtr_;

        // Not necessary to resize the list as there are no retired points
        // directAddr.setSize(pMesh_.size());

        insertedPointLabelsPtr_ = new labelList(pMesh_.size());
        labelList& insertedPoints = *insertedPointLabelsPtr_;

        label nInsertedPoints = 0;

        forAll(directAddr, pointI)
        {
            if (directAddr[pointI] < 0)
            {
                // Found inserted point
                directAddr[pointI] = 0;
                insertedPoints[nInsertedPoints] = pointI;
                nInsertedPoints++;
            }
        }

        insertedPoints.setSize(nInsertedPoints);
    }
    else
    {
        // Interpolative addressing

        interpolationAddrPtr_ = new labelListList(pMesh_.size());
        labelListList& addr = *interpolationAddrPtr_;

        weightsPtr_ = new scalarListList(pMesh_.size());
        scalarListList& w = *weightsPtr_;

        // Points created from other points (i.e. points merged into it).
        const List<objectMap>& cfc = mpm_.pointsFromPointsMap();

        forAll(cfc, cfcI)
        {
            // Get addressing
            const labelList& mo = cfc[cfcI].masterObjects();

            label pointI = cfc[cfcI].index();

            if (addr[pointI].size())
            {
                FatalErrorIn("void pointMapper::calcAddressing() const")
                    << "Master point " << pointI
                    << " mapped from points " << mo
                    << " already destination of mapping." << abort(FatalError);
            }

            // Map from masters, uniform weights
            addr[pointI] = mo;
            w[pointI] = scalarList(mo.size(), 1.0/mo.size());
        }


        // Do mapped points. Note that can already be set from pointsFromPoints
        // so check if addressing size still zero.

        const labelList& cm = mpm_.pointMap();

        forAll(cm, pointI)
        {
            if (cm[pointI] > -1 && addr[pointI].empty())
            {
                // Mapped from a single point
                addr[pointI] = labelList(1, cm[pointI]);
                w[pointI] = scalarList(1, 1.0);
            }
        }

        // Grab inserted points (for them the size of addressing is still zero)

        insertedPointLabelsPtr_ = new labelList(pMesh_.size());
        labelList& insertedPoints = *insertedPointLabelsPtr_;

        label nInsertedPoints = 0;

        forAll(addr, pointI)
        {
            if (addr[pointI].empty())
            {
                // Mapped from a dummy point. Take point 0 with weight 1.
                addr[pointI] = labelList(1, label(0));
                w[pointI] = scalarList(1, 1.0);

                insertedPoints[nInsertedPoints] = pointI;
                nInsertedPoints++;
            }
        }

        insertedPoints.setSize(nInsertedPoints);
    }
}


void Foam::pointMapper::clearOut()
{
    deleteDemandDrivenData(directAddrPtr_);
    deleteDemandDrivenData(interpolationAddrPtr_);
    deleteDemandDrivenData(weightsPtr_);
    deleteDemandDrivenData(insertedPointLabelsPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::pointMapper::pointMapper(const pointMesh& pMesh, const mapPolyMesh& mpm)
:
    pMesh_(pMesh),
    mpm_(mpm),
    insertedPoints_(true),
    direct_(false),
    directAddrPtr_(NULL),
    interpolationAddrPtr_(NULL),
    weightsPtr_(NULL),
    insertedPointLabelsPtr_(NULL)
{
    // Check for possibility of direct mapping
    if (mpm_.pointsFromPointsMap().empty())
    {
        direct_ = true;
    }
    else
    {
        direct_ = false;
    }

    // Check for inserted points
    if (direct_ && (mpm_.pointMap().empty() || min(mpm_.pointMap()) > -1))
    {
        insertedPoints_ = false;
    }
    else
    {
        //Check if there are inserted points with no owner

        // Make a copy of the point map, add the entries for points from points
        // and check for left-overs
        labelList cm(pMesh_.size(), -1);

        const List<objectMap>& cfc = mpm_.pointsFromPointsMap();

        forAll(cfc, cfcI)
        {
            cm[cfc[cfcI].index()] = 0;
        }

        if (min(cm) < 0)
        {
            insertedPoints_ = true;
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pointMapper::~pointMapper()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::pointMapper::size() const
{
    return mpm_.pointMap().size();
}


Foam::label Foam::pointMapper::sizeBeforeMapping() const
{
    return mpm_.nOldPoints();
}


const Foam::labelUList& Foam::pointMapper::directAddressing() const
{
    if (!direct())
    {
        FatalErrorIn
        (
            "const labelUList& pointMapper::directAddressing() const"
        )   << "Requested direct addressing for an interpolative mapper."
            << abort(FatalError);
    }

    if (!insertedObjects())
    {
        // No inserted points.  Re-use pointMap
        return mpm_.pointMap();
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


const Foam::labelListList& Foam::pointMapper::addressing() const
{
    if (direct())
    {
        FatalErrorIn
        (
            "const labelListList& pointMapper::addressing() const"
        )   << "Requested interpolative addressing for a direct mapper."
            << abort(FatalError);
    }

    if (!interpolationAddrPtr_)
    {
        calcAddressing();
    }

    return *interpolationAddrPtr_;
}


const Foam::scalarListList& Foam::pointMapper::weights() const
{
    if (direct())
    {
        FatalErrorIn
        (
            "const scalarListList& pointMapper::weights() const"
        )   << "Requested interpolative weights for a direct mapper."
            << abort(FatalError);
    }

    if (!weightsPtr_)
    {
        calcAddressing();
    }

    return *weightsPtr_;
}


const Foam::labelList& Foam::pointMapper::insertedObjectLabels() const
{
    if (!insertedPointLabelsPtr_)
    {
        if (!insertedObjects())
        {
            // There are no inserted points
            insertedPointLabelsPtr_ = new labelList(0);
        }
        else
        {
            calcAddressing();
        }
    }

    return *insertedPointLabelsPtr_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// ************************************************************************* //
