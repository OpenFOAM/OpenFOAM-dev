/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellEdgeAddressing::cellEdgeAddressing
(
    const cell& c,
    const faceList& fs,
    const bool cOwnsFirst
)
{
    // Compute the number of cell edges
    label nCellEdges = 0;
    forAll(c, cfi)
    {
        const face& f = fs[c[cfi]];

        nCellEdges += f.size();
    }
    nCellEdges /= 2;

    // Construct a map to enumerate the cell edges
    HashList<label, edge, QuickHashEdge> edgeToCei(nCellEdges*6);
    {
        label cellEdgei = 0;

        forAll(c, cfi)
        {
            const face& f = fs[c[cfi]];

            forAll(f, fei)
            {
                cellEdgei += edgeToCei.insert(f.faceEdge(fei), cellEdgei);
            }
        }
    }

    // Allocate and initialise the addressing
    cfiAndFeiToCei_ = labelListList(c.size());
    ceiToCfiAndFei_ =
        List<Pair<labelPair>>
        (
            nCellEdges,
            Pair<labelPair>(labelPair(-1, -1), labelPair(-1, -1))
        );

    // Copy data out of the map into the addressing
    forAll(c, cfi)
    {
        const face& f = fs[c[cfi]];

        cfiAndFeiToCei_[cfi] = labelList(f.size(), -1);

        forAll(f, fei)
        {
            const label cei = edgeToCei[f.faceEdge(fei)];

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
    boolList cOwnsIsSet(c.size(), false);
    cOwnsIsSet[0] = true;

    // Compare cell-face-edges to determine face signs
    forAll(c, cfi)
    {
        const face& f = fs[c[cfi]];

        if (!cOwnsIsSet[cfi]) continue;

        forAll(f, fei)
        {
            const label cei = cfiAndFeiToCei_[cfi][fei];

            const labelPair& other =
                ceiToCfiAndFei_[cei]
                [
                    ceiToCfiAndFei_[cei][0] == labelPair(cfi, fei)
                ];
            const label cfj = other[0], fej = other[1];

            if (cOwnsIsSet[cfj]) continue;

            const label sign =
                edge::compare
                (
                    f.faceEdge(fei),
                    fs[c[cfj]].faceEdge(fej)
                );

            cOwns_[cfj] = sign < 0 ? cOwns_[cfi] : !cOwns_[cfi];
            cOwnsIsSet[cfj] = true;
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


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

const Foam::cellEdgeAddressing& Foam::cellEdgeAddressingList::operator[]
(
    const label celli
) const
{
    if (!list_.set(celli))
    {
        const cell& c = mesh().cells()[celli];
        const faceList& fs = mesh().faces();
        const labelList& fOwners = mesh().faceOwner();

        list_.set
        (
            celli,
            new cellEdgeAddressing(c, fs, fOwners[c[0]] == celli)
        );
    }

    return list_[celli];
}


// ************************************************************************* //
