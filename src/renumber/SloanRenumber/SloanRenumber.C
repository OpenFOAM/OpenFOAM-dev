/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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

    Adapted from Boost graph/example/sloan_ordering.cpp

\*---------------------------------------------------------------------------*/

#include "SloanRenumber.H"
#include "addToRunTimeSelectionTable.H"
#include "decompositionMethod.H"
#include "processorPolyPatch.H"
#include "syncTools.H"

#include <boost/config.hpp>
#include <vector>
#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/sloan_ordering.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/bandwidth.hpp>
#include <boost/graph/profile.hpp>
#include <boost/graph/wavefront.hpp>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace boost;
using namespace std;

//Defining the graph type
typedef adjacency_list
<
    setS,
    vecS,
    undirectedS,
    property
    <
        vertex_color_t,
        default_color_type,
        property
        <
            vertex_degree_t,
            Foam::label,
            property
            <
                vertex_priority_t,
                Foam::scalar
            >
        >
    >
> Graph;

typedef graph_traits<Graph>::vertex_descriptor Vertex;
typedef graph_traits<Graph>::vertices_size_type size_type;

namespace Foam
{
    defineTypeNameAndDebug(SloanRenumber, 0);

    addToRunTimeSelectionTable
    (
        renumberMethod,
        SloanRenumber,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SloanRenumber::SloanRenumber(const dictionary& renumberDict)
:
    renumberMethod(renumberDict),
    reverse_
    (
        renumberDict.optionalSubDict
        (
            typeName + "Coeffs"
        ).lookupOrDefault<Switch>("reverse", false)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::SloanRenumber::renumber
(
    const polyMesh& mesh,
    const pointField& points
) const
{
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    // Construct graph : faceOwner + connections across cyclics.

    // Determine neighbour cell
    labelList nbr(mesh.nFaces()-mesh.nInternalFaces(), -1);
    forAll(pbm, patchi)
    {
        if (pbm[patchi].coupled() && !isA<processorPolyPatch>(pbm[patchi]))
        {
            SubList<label>
            (
                nbr,
                pbm[patchi].size(),
                pbm[patchi].start()-mesh.nInternalFaces()
            ) = pbm[patchi].faceCells();
        }
    }
    syncTools::swapBoundaryFaceList(mesh, nbr);


    Graph G(mesh.nCells());

    // Add internal faces
    forAll(mesh.faceNeighbour(), facei)
    {
        add_edge(mesh.faceOwner()[facei], mesh.faceNeighbour()[facei], G);
    }
    // Add cyclics
    forAll(pbm, patchi)
    {
        if
        (
            pbm[patchi].coupled()
        && !isA<processorPolyPatch>(pbm[patchi])
        &&  refCast<const coupledPolyPatch>(pbm[patchi]).owner()
        )
        {
            const labelUList& faceCells = pbm[patchi].faceCells();
            forAll(faceCells, i)
            {
                label bFacei = pbm[patchi].start()+i-mesh.nInternalFaces();
                label nbrCelli = nbr[bFacei];

                if (faceCells[i] < nbrCelli)
                {
                    add_edge(faceCells[i], nbrCelli, G);
                }
                else
                {
                    add_edge(nbrCelli, faceCells[i], G);
                }
            }
        }
    }


    // Creating two iterators over the vertices
    graph_traits<Graph>::vertex_iterator ui, ui_end;

    // Creating a property_map with the degrees of the degrees of each vertex
    property_map<Graph,vertex_degree_t>::type deg = get(vertex_degree, G);
    for (boost::tie(ui, ui_end) = vertices(G); ui != ui_end; ++ui)
        deg[*ui] = degree(*ui, G);

    // Creating a property_map for the indices of a vertex
    property_map<Graph, vertex_index_t>::type index_map = get(vertex_index, G);

    // Creating a vector of vertices
    std::vector<Vertex> sloan_order(num_vertices(G));

    sloan_ordering
    (
        G,
        sloan_order.begin(),
        get(vertex_color, G),
        make_degree_map(G),
        get(vertex_priority, G)
    );

    labelList orderedToOld(sloan_order.size());
    forAll(orderedToOld, c)
    {
        orderedToOld[c] = index_map[sloan_order[c]];
    }

    if (reverse_)
    {
        reverse(orderedToOld);
    }

    return orderedToOld;
}


Foam::labelList Foam::SloanRenumber::renumber
(
    const labelListList& cellCells,
    const pointField& points
) const
{
    Graph G(cellCells.size());

    forAll(cellCells, celli)
    {
        const labelList& nbrs = cellCells[celli];
        forAll(nbrs, i)
        {
            if (nbrs[i] > celli)
            {
                add_edge(celli, nbrs[i], G);
            }
        }
    }

    // Creating two iterators over the vertices
    graph_traits<Graph>::vertex_iterator ui, ui_end;

    // Creating a property_map with the degrees of the degrees of each vertex
    property_map<Graph,vertex_degree_t>::type deg = get(vertex_degree, G);
    for (boost::tie(ui, ui_end) = vertices(G); ui != ui_end; ++ui)
        deg[*ui] = degree(*ui, G);

    // Creating a property_map for the indices of a vertex
    property_map<Graph, vertex_index_t>::type index_map = get(vertex_index, G);

    // Creating a vector of vertices
    std::vector<Vertex> sloan_order(num_vertices(G));

    sloan_ordering
    (
        G,
        sloan_order.begin(),
        get(vertex_color, G),
        make_degree_map(G),
        get(vertex_priority, G)
    );

    labelList orderedToOld(sloan_order.size());
    forAll(orderedToOld, c)
    {
        orderedToOld[c] = index_map[sloan_order[c]];
    }

    if (reverse_)
    {
        reverse(orderedToOld);
    }

    return orderedToOld;
}


// ************************************************************************* //
