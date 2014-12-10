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

#include "primitiveMesh.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelListList& Foam::primitiveMesh::cellPoints() const
{
    if (!cpPtr_)
    {
        if (debug)
        {
            Pout<< "primitiveMesh::cellPoints() : "
                << "calculating cellPoints" << endl;

            if (debug == -1)
            {
                // For checking calls:abort so we can quickly hunt down
                // origin of call
                FatalErrorIn("primitiveMesh::cellPoints()")
                    << abort(FatalError);
            }
        }

        // Invert pointCells
        cpPtr_ = new labelListList(nCells());
        invertManyToMany(nCells(), pointCells(), *cpPtr_);
    }

    return *cpPtr_;
}


const Foam::labelList& Foam::primitiveMesh::cellPoints
(
    const label cellI,
    DynamicList<label>& storage
) const
{
    if (hasCellPoints())
    {
        return cellPoints()[cellI];
    }
    else
    {
        const faceList& fcs = faces();
        const labelList& cFaces = cells()[cellI];

        labelSet_.clear();

        forAll(cFaces, i)
        {
            const labelList& f = fcs[cFaces[i]];

            forAll(f, fp)
            {
                labelSet_.insert(f[fp]);
            }
        }

        storage.clear();
        if (labelSet_.size() > storage.capacity())
        {
            storage.setCapacity(labelSet_.size());
        }

        forAllConstIter(labelHashSet, labelSet_, iter)
        {
            storage.append(iter.key());
        }

        return storage;
    }
}


const Foam::labelList& Foam::primitiveMesh::cellPoints(const label cellI) const
{
    return cellPoints(cellI, labels_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
