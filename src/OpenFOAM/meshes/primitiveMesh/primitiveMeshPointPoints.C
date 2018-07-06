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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::primitiveMesh::calcPointPoints() const
{
    if (debug)
    {
        Pout<< "primitiveMesh::calcPointPoints() : "
            << "calculating pointPoints"
            << endl;

        if (debug == -1)
        {
            // For checking calls:abort so we can quickly hunt down
            // origin of call
            FatalErrorInFunction
                << abort(FatalError);
        }
    }

    // It is an error to attempt to recalculate pointPoints
    // if the pointer is already set
    if (ppPtr_)
    {
        FatalErrorInFunction
            << "pointPoints already calculated"
            << abort(FatalError);
    }
    else
    {
        const edgeList& e = edges();
        const labelListList& pe = pointEdges();

        ppPtr_ = new labelListList(pe.size());
        labelListList& pp = *ppPtr_;

        forAll(pe, pointi)
        {
            pp[pointi].setSize(pe[pointi].size());

            forAll(pe[pointi], ppi)
            {
                if (e[pe[pointi][ppi]].start() == pointi)
                {
                    pp[pointi][ppi] = e[pe[pointi][ppi]].end();
                }
                else if (e[pe[pointi][ppi]].end() == pointi)
                {
                    pp[pointi][ppi] = e[pe[pointi][ppi]].start();
                }
                else
                {
                    FatalErrorInFunction
                        << "something wrong with edges"
                        << abort(FatalError);
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelListList& Foam::primitiveMesh::pointPoints() const
{
    if (!ppPtr_)
    {
        calcPointPoints();
    }

    return *ppPtr_;
}


const Foam::labelList& Foam::primitiveMesh::pointPoints
(
    const label pointi,
    DynamicList<label>& storage
) const
{
    if (hasPointPoints())
    {
        return pointPoints()[pointi];
    }
    else
    {
        const edgeList& edges = this->edges();
        const labelList& pEdges = pointEdges()[pointi];

        storage.clear();

        if (pEdges.size() > storage.capacity())
        {
            storage.setCapacity(pEdges.size());
        }

        forAll(pEdges, i)
        {
            storage.append(edges[pEdges[i]].otherVertex(pointi));
        }

        return storage;
    }
}


const Foam::labelList& Foam::primitiveMesh::pointPoints
(
    const label pointi
) const
{
    return pointPoints(pointi, labels_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
