/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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

#include "primitiveMesh.H"
#include "cell.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::primitiveMesh::calcPointCells() const
{
    // Loop through cells and mark up points

    if (debug)
    {
        Pout<< "primitiveMesh::calcPointCells() : "
            << "calculating pointCells"
            << endl;

        if (debug == -1)
        {
            // For checking calls:abort so we can quickly hunt down
            // origin of call
            FatalErrorInFunction
                << abort(FatalError);
        }
    }

    // It is an error to attempt to recalculate pointCells
    // if the pointer is already set
    if (pcPtr_)
    {
        FatalErrorInFunction
            << "pointCells already calculated"
            << abort(FatalError);
    }
    else
    {
        const cellList& cf = cells();

        // Count number of cells per point

        labelList npc(nPoints(), 0);

        forAll(cf, celli)
        {
            const labelList curPoints = cf[celli].labels(faces());

            forAll(curPoints, pointi)
            {
                label ptI = curPoints[pointi];

                npc[ptI]++;
            }
        }


        // Size and fill cells per point

        pcPtr_ = new labelListList(npc.size());
        labelListList& pointCellAddr = *pcPtr_;

        forAll(pointCellAddr, pointi)
        {
            pointCellAddr[pointi].setSize(npc[pointi]);
        }
        npc = 0;


        forAll(cf, celli)
        {
            const labelList curPoints = cf[celli].labels(faces());

            forAll(curPoints, pointi)
            {
                label ptI = curPoints[pointi];

                pointCellAddr[ptI][npc[ptI]++] = celli;
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelListList& Foam::primitiveMesh::pointCells() const
{
    if (!pcPtr_)
    {
        calcPointCells();
    }

    return *pcPtr_;
}


const Foam::labelList& Foam::primitiveMesh::pointCells
(
    const label pointi,
    DynamicList<label>& storage
) const
{
    if (hasPointCells())
    {
        return pointCells()[pointi];
    }
    else
    {
        const labelList& own = faceOwner();
        const labelList& nei = faceNeighbour();
        const labelList& pFaces = pointFaces()[pointi];

        storage.clear();

        forAll(pFaces, i)
        {
            const label facei = pFaces[i];

            // Append owner
            storage.append(own[facei]);

            // Append neighbour
            if (facei < nInternalFaces())
            {
                storage.append(nei[facei]);
            }
        }

        // Filter duplicates
        if (storage.size() > 1)
        {
            sort(storage);

            label n = 1;
            for (label i = 1; i < storage.size(); i++)
            {
                if (storage[i-1] != storage[i])
                {
                    storage[n++] = storage[i];
                }
            }

            // truncate addressed list
            storage.setSize(n);
        }

        return storage;
    }
}


const Foam::labelList& Foam::primitiveMesh::pointCells(const label pointi) const
{
    return pointCells(pointi, labels_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
