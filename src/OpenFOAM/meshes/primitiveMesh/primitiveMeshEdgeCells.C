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
#include "ListOps.H"


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelListList& Foam::primitiveMesh::edgeCells() const
{
    if (!ecPtr_)
    {
        if (debug)
        {
            Pout<< "primitiveMesh::edgeCells() : calculating edgeCells" << endl;

            if (debug == -1)
            {
                // For checking calls:abort so we can quickly hunt down
                // origin of call
                FatalErrorInFunction
                    << abort(FatalError);
            }
        }
        // Invert cellEdges
        ecPtr_ = new labelListList(nEdges());
        invertManyToMany(nEdges(), cellEdges(), *ecPtr_);
    }

    return *ecPtr_;
}


const Foam::labelList& Foam::primitiveMesh::edgeCells
(
    const label edgeI,
    DynamicList<label>& storage
) const
{
    if (hasEdgeCells())
    {
        return edgeCells()[edgeI];
    }
    else
    {
        const labelList& own = faceOwner();
        const labelList& nei = faceNeighbour();

        // Construct edgeFaces
        DynamicList<label> eFacesStorage;
        const labelList& eFaces = edgeFaces(edgeI, eFacesStorage);

        storage.clear();

        // Do quadratic insertion.
        forAll(eFaces, i)
        {
            label facei = eFaces[i];

            {
                label ownCelli = own[facei];

                // Check if not already in storage
                forAll(storage, j)
                {
                    if (storage[j] == ownCelli)
                    {
                        ownCelli = -1;
                        break;
                    }
                }

                if (ownCelli != -1)
                {
                    storage.append(ownCelli);
                }
            }

            if (isInternalFace(facei))
            {
                label neiCelli = nei[facei];

                forAll(storage, j)
                {
                    if (storage[j] == neiCelli)
                    {
                        neiCelli = -1;
                        break;
                    }
                }

                if (neiCelli != -1)
                {
                    storage.append(neiCelli);
                }
            }
        }

        return storage;
    }
}


const Foam::labelList& Foam::primitiveMesh::edgeCells(const label edgeI) const
{
    return edgeCells(edgeI, labels_);
}


// ************************************************************************* //
