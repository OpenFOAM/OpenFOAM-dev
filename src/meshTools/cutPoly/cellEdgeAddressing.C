/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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

#include "cellEdgeAddressing.H"
#include "polyDistributionMap.H"
#include "polyTopoChangeMap.H"
#include "polyMeshMap.H"
#include "HashList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// An edge hash functor that is much faster than Hash<edge>, but is also more
// prone to collisions
struct QuickHashEdge
{
    unsigned operator()(const edge& e)
    {
        return e[0]*e[1];
    }
};

// An invalid null edge for unset entries in the edge map
template<>
const edge HashList<label, edge, QuickHashEdge>::nullKey(-labelMax, -labelMax);

// Workspace for addressing calculations
struct cellEdgeAddressingWorkspace
{
    HashList<label, edge, QuickHashEdge> edgeToCei;

    DynamicList<label> cfei0;

    cellEdgeAddressingWorkspace()
    :
        edgeToCei(0),
        cfei0()
    {
        resizeAndClear(6, 12);
    }

    void resizeAndClear(const label nCellFaces, const label nCellEdges)
    {
        if (edgeToCei.capacity() < nCellEdges)
        {
            edgeToCei.resizeAndClear(nCellEdges*6);
        }
        else
        {
            edgeToCei.clear();
        }

        cfei0.resize(nCellFaces);
        cfei0 = -1;
    }
}
workspace_;

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellEdgeAddressingData::cellEdgeAddressingData
(
    const cell& c,
    const faceList& fs,
    const bool cOwnsFirst
)
{
    // Allocate and initialise the addressing
    cfiAndFeiToCei_.resize(UIndirectList<face>(fs, c));
    ceiToCfiAndFei_ =
        List<Pair<labelPair>>
        (
            cfiAndFeiToCei_.m().size()/2,
            Pair<labelPair>(labelPair(-1, -1), labelPair(-1, -1))
        );

    // Resize and reset the workspace
    workspace_.resizeAndClear(c.size(), ceiToCfiAndFei_.size());

    // Construct a map to enumerate the cell edges
    {
        label cellEdgei = 0;

        forAll(c, cfi)
        {
            const face& f = fs[c[cfi]];

            forAll(f, fei)
            {
                cellEdgei +=
                    workspace_.edgeToCei.insert(f.faceEdge(fei), cellEdgei);
            }
        }
    }

    // Copy data out of the map into the addressing
    forAll(c, cfi)
    {
        const face& f = fs[c[cfi]];

        cfiAndFeiToCei_[cfi] = -1;

        forAll(f, fei)
        {
            const label cei = workspace_.edgeToCei[f.faceEdge(fei)];

            cfiAndFeiToCei_[cfi][fei] = cei;
            ceiToCfiAndFei_[cei]
            [
                ceiToCfiAndFei_[cei][0] != labelPair(-1, -1)
            ] = {cfi, fei};
        }
    }

    // Allocate and initialise the face ownership
    cOwns_ = boolList(c.size(), false);
    cOwns_[0] = cOwnsFirst;
    workspace_.cfei0[0] = 0;

    // Walk around the cell, comparing edges to determine face ownership
    {
        label cfi = 0, fei = 0;

        do
        {
            // Get the adjacent face and face-edge (given subscript j)
            const label cei = cfiAndFeiToCei_[cfi][fei];
            const labelPair cfiAndFei(labelPair(cfi, fei));
            const labelPair& cfjAndFej =
                ceiToCfiAndFei_[cei][ceiToCfiAndFei_[cei][0] == cfiAndFei];
            const label cfj = cfjAndFej[0], fej = cfjAndFej[1];

            // If the adjacent face has not been visited then set its ownership
            // and it's starting face edge and move forwards into it
            if (workspace_.cfei0[cfj] == -1)
            {
                // If the face-edges point in different directions then the
                // faces have the same owner. If they point in the same
                // direction then they have different owners.
                const label sign =
                    edge::compare
                    (
                        fs[c[cfi]].faceEdge(fei),
                        fs[c[cfj]].faceEdge(fej)
                    );

                cOwns_[cfj] = sign < 0 ? cOwns_[cfi] : !cOwns_[cfi];

                workspace_.cfei0[cfj] = fej;
                cfi = cfj;
                fei = (fej + 1) % fs[c[cfj]].size();
            }

            // If the adjacent face has been visited, and we are not back at
            // the starting edge of this face, then move to the next edge of
            // this face
            else if (fei != workspace_.cfei0[cfi])
            {
                fei = (fei + 1) % fs[c[cfi]].size();
            }

            // If we are back at the first edge then this face is complete, so
            // move backwards into the adjacent face
            else // if (fei == workspace_.cfei0[cfi])
            {
                cfi = cfj;
                fei = fej;
            }
        }
        while (cfi != 0 || fei != 0);
    }
}


Foam::cellEdgeAddressingList::cellEdgeAddressingList(const polyMesh& mesh)
:
    DemandDrivenMeshObject
    <
        polyMesh,
        UpdateableMeshObject,
        cellEdgeAddressingList
    >(mesh),
    list_(mesh.nCells())
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::cellEdgeAddressingList::movePoints()
{
    return true;
}


void Foam::cellEdgeAddressingList::distribute(const polyDistributionMap& map)
{
    // We could be more selective here and try to keep addressing for cells
    // that haven't changed

    list_.clear();
    list_.resize(map.mesh().nCells());
}


void Foam::cellEdgeAddressingList::topoChange(const polyTopoChangeMap& map)
{
    // We could be more selective here and try to keep addressing for cells
    // that haven't changed

    list_.clear();
    list_.resize(map.mesh().nCells());
}


void Foam::cellEdgeAddressingList::mapMesh(const polyMeshMap& map)
{
    list_.clear();
    list_.resize(map.mesh().nCells());
}


// ************************************************************************* //
