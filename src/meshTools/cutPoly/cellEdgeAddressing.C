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

    DynamicList<bool> cOwnsIsSet;

    cellEdgeAddressingWorkspace()
    :
        edgeToCei(0),
        cOwnsIsSet()
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

        cOwnsIsSet.resize(nCellFaces);
        cOwnsIsSet = false;
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

    // Allocate and initialise the face signs
    cOwns_ = boolList(c.size(), false);
    cOwns_[0] = cOwnsFirst;
    workspace_.cOwnsIsSet[0] = true;

    // Compare cell-face-edges to determine face signs
    forAll(c, cfi)
    {
        const face& f = fs[c[cfi]];

        if (!workspace_.cOwnsIsSet[cfi]) continue;

        forAll(f, fei)
        {
            const label cei = cfiAndFeiToCei_[cfi][fei];

            const labelPair& other =
                ceiToCfiAndFei_[cei]
                [
                    ceiToCfiAndFei_[cei][0] == labelPair(cfi, fei)
                ];
            const label cfj = other[0], fej = other[1];

            if (workspace_.cOwnsIsSet[cfj]) continue;

            const label sign =
                edge::compare
                (
                    f.faceEdge(fei),
                    fs[c[cfj]].faceEdge(fej)
                );

            cOwns_[cfj] = sign < 0 ? cOwns_[cfi] : !cOwns_[cfi];
            workspace_.cOwnsIsSet[cfj] = true;
        }
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
