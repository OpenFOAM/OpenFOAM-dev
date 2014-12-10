/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

Class
    Foam::syncTools

Description
    Various tools to aid synchronizing lists across coupled patches. WIP.

    Require
    - combineOperator (e.g. sumEqOp - not sumOp!) that is defined for the
      type and combineReduce(UList\<T\>, combineOperator) should be defined.
    - null value which gets overridden by any valid value.
    - transform function

SourceFiles
    syncTools.C
    syncToolsTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef syncTools_H
#define syncTools_H

#include "Pstream.H"
#include "EdgeMap.H"
#include "PackedBoolList.H"
#include "polyMesh.H"
#include "coupledPolyPatch.H"
#include "mapDistribute.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyBoundaryMesh;

/*---------------------------------------------------------------------------*\
                           Class syncTools Declaration
\*---------------------------------------------------------------------------*/

class syncTools
{
    // Private Member Functions

        //- Combine value with existing value in map.
        template<class T, class CombineOp>
        static void combine
        (
            Map<T>& pointValues,
            const CombineOp& cop,
            const label index,
            const T& val
        );

        //- Combine val with existing value at (implicit index) e.
        template<class T, class CombineOp>
        static void combine
        (
            EdgeMap<T>& edgeValues,
            const CombineOp& cop,
            const edge& index,
            const T& val
        );


public:


        // Basic routines with user-supplied transformation. Preferably
        // use specialisations below.

            //- Synchronize values on selected points.
            template<class T, class CombineOp, class TransformOp>
            static void syncPointMap
            (
                const polyMesh&,
                Map<T>& pointValues,
                const CombineOp& cop,
                const TransformOp& top
            );

            //- Synchronize values on selected edges.
            template<class T, class CombineOp, class TransformOp>
            static void syncEdgeMap
            (
                const polyMesh&,
                EdgeMap<T>& edgeValues,
                const CombineOp& cop,
                const TransformOp& top
            );

            //- Synchronize values on all mesh points.
            template<class T, class CombineOp, class TransformOp>
            static void syncPointList
            (
                const polyMesh&,
                List<T>&,
                const CombineOp& cop,
                const T& nullValue,
                const TransformOp& top
            );

            //- Synchronize values on selected mesh points.
            template<class T, class CombineOp, class TransformOp>
            static void syncPointList
            (
                const polyMesh&,
                const labelList& meshPoints,
                List<T>&,
                const CombineOp& cop,
                const T& nullValue,
                const TransformOp& top
            );

            //- Synchronize values on all mesh edges.
            template<class T, class CombineOp, class TransformOp>
            static void syncEdgeList
            (
                const polyMesh&,
                List<T>&,
                const CombineOp& cop,
                const T& nullValue,
                const TransformOp& top
            );

            //- Synchronize values on selected mesh edges.
            template<class T, class CombineOp, class TransformOp>
            static void syncEdgeList
            (
                const polyMesh&,
                const labelList& meshEdges,
                List<T>&,
                const CombineOp& cop,
                const T& nullValue,
                const TransformOp& top
            );

            //- Synchronize values on boundary faces only.
            template<class T, class CombineOp, class TransformOp>
            static void syncBoundaryFaceList
            (
                const polyMesh&,
                UList<T>&,
                const CombineOp& cop,
                const TransformOp& top
            );


        // Synchronise point-wise data

            //- Synchronize values on all mesh points.
            template<class T, class CombineOp>
            static void syncPointList
            (
                const polyMesh& mesh,
                List<T>& l,
                const CombineOp& cop,
                const T& nullValue
            )
            {
                syncPointList
                (
                    mesh,
                    l,
                    cop,
                    nullValue,
                    mapDistribute::transform()
                );
            }

            //- Synchronize locations on all mesh points.
            template<class CombineOp>
            static void syncPointPositions
            (
                const polyMesh& mesh,
                List<point>& l,
                const CombineOp& cop,
                const point& nullValue
            )
            {
                syncPointList
                (
                    mesh,
                    l,
                    cop,
                    nullValue,
                    mapDistribute::transformPosition()
                );
            }

            //- Synchronize values on selected mesh points.
            template<class T, class CombineOp>
            static void syncPointList
            (
                const polyMesh& mesh,
                const labelList& meshPoints,
                List<T>& l,
                const CombineOp& cop,
                const T& nullValue
            )
            {
                syncPointList
                (
                    mesh,
                    meshPoints,
                    l,
                    cop,
                    nullValue,
                    mapDistribute::transform()
                );
            }

            //- Synchronize locations on selected mesh points.
            template<class CombineOp>
            static void syncPointPositions
            (
                const polyMesh& mesh,
                const labelList& meshPoints,
                List<point>& l,
                const CombineOp& cop,
                const point& nullValue
            )
            {
                syncPointList
                (
                    mesh,
                    meshPoints,
                    l,
                    cop,
                    nullValue,
                    mapDistribute::transformPosition()
                );
            }


        // Synchronise edge-wise data

            //- Synchronize values on all mesh edges.
            template<class T, class CombineOp>
            static void syncEdgeList
            (
                const polyMesh& mesh,
                List<T>& l,
                const CombineOp& cop,
                const T& nullValue
            )
            {
                syncEdgeList
                (
                    mesh,
                    l,
                    cop,
                    nullValue,
                    mapDistribute::transform()
                );
            }

            //- Synchronize locations on all mesh edges.
            template<class CombineOp>
            static void syncEdgePositions
            (
                const polyMesh& mesh,
                List<point>& l,
                const CombineOp& cop,
                const point& nullValue
            )
            {
                syncEdgeList
                (
                    mesh,
                    l,
                    cop,
                    nullValue,
                    mapDistribute::transformPosition()
                );
            }

            //- Synchronize values on selected mesh edges.
            template<class T, class CombineOp>
            static void syncEdgeList
            (
                const polyMesh& mesh,
                const labelList& meshEdges,
                List<T>& l,
                const CombineOp& cop,
                const T& nullValue
            )
            {
                syncEdgeList
                (
                    mesh,
                    meshEdges,
                    l,
                    cop,
                    nullValue,
                    mapDistribute::transform()
                );
            }

            //- Synchronize locations on selected mesh edges.
            template<class CombineOp>
            static void syncEdgePositions
            (
                const polyMesh& mesh,
                const labelList& meshEdges,
                List<point>& l,
                const CombineOp& cop,
                const point& nullValue
            )
            {
                syncEdgeList
                (
                    mesh,
                    meshEdges,
                    l,
                    cop,
                    nullValue,
                    mapDistribute::transformPosition()
                );
            }



        // Synchronise face-wise data

            //- Synchronize values on boundary faces only.
            template<class T, class CombineOp>
            static void syncBoundaryFaceList
            (
                const polyMesh& mesh,
                UList<T>& l,
                const CombineOp& cop
            )
            {
                syncBoundaryFaceList(mesh, l, cop, mapDistribute::transform());
            }

            //- Synchronize locations on boundary faces only.
            template<class CombineOp>
            static void syncBoundaryFacePositions
            (
                const polyMesh& mesh,
                UList<point>& l,
                const CombineOp& cop
            )
            {
                syncBoundaryFaceList
                (
                    mesh,
                    l,
                    cop,
                    mapDistribute::transformPosition()
                );
            }

            //- Synchronize values on all mesh faces.
            template<class T, class CombineOp>
            static void syncFaceList
            (
                const polyMesh& mesh,
                UList<T>& l,
                const CombineOp& cop
            )
            {
                SubList<T> bndValues
                (
                    l,
                    mesh.nFaces()-mesh.nInternalFaces(),
                    mesh.nInternalFaces()
                );

                syncBoundaryFaceList
                (
                    mesh,
                    bndValues,
                    cop,
                    mapDistribute::transform()
                );
            }

            //- Synchronize locations on all mesh faces.
            template<class CombineOp>
            static void syncFacePositions
            (
                const polyMesh& mesh,
                UList<point>& l,
                const CombineOp& cop
            )
            {
                SubList<point> bndValues
                (
                    l,
                    mesh.nFaces()-mesh.nInternalFaces(),
                    mesh.nInternalFaces()
                );
                syncBoundaryFaceList
                (
                    mesh,
                    bndValues,
                    cop,
                    mapDistribute::transformPosition()
                );
            }

            //- Swap coupled boundary face values.
            template<class T>
            static void swapBoundaryFaceList
            (
                const polyMesh& mesh,
                UList<T>& l
            )
            {
                syncBoundaryFaceList
                (
                    mesh,
                    l,
                    eqOp<T>(),
                    mapDistribute::transform()
                );
            }

             //- Swap coupled positions.
            static void swapBoundaryFacePositions
            (
                const polyMesh& mesh,
                UList<point>& l
            )
            {
                syncBoundaryFaceList
                (
                    mesh,
                    l,
                    eqOp<point>(),
                    mapDistribute::transformPosition()
                );
            }

            //- Swap coupled face values.
            template<class T>
            static void swapFaceList
            (
                const polyMesh& mesh,
                UList<T>& l
            )
            {
                SubList<T> bndValues
                (
                    l,
                    mesh.nFaces()-mesh.nInternalFaces(),
                    mesh.nInternalFaces()
                );
                syncBoundaryFaceList
                (
                    mesh,
                    bndValues,
                    eqOp<T>(),
                    mapDistribute::transform()
                );
            }

            //- Swap to obtain neighbour cell values for all boundary faces
            template<class T>
            static void swapBoundaryCellList
            (
                const polyMesh& mesh,
                const UList<T>& cellData,
                List<T>& neighbourCellData
            );

            //- Swap to obtain neighbour cell positions for all boundary faces
            static void swapBoundaryCellPositions
            (
                const polyMesh& mesh,
                const UList<point>& cellData,
                List<point>& neighbourCellData
            );

        // Sparse versions

            //- Synchronize values on selected points.
            template<class T, class CombineOp>
            static void syncPointMap
            (
                const polyMesh& mesh,
                Map<T>& l,
                const CombineOp& cop
            )
            {
                syncPointMap(mesh, l, cop, mapDistribute::transform());
            }

            //- Synchronize locations on selected points.
            template<class CombineOp>
            static void syncPointPositions
            (
                const polyMesh& mesh,
                Map<point>& l,
                const CombineOp& cop
            )
            {
                syncPointMap(mesh, l, cop, mapDistribute::transformPosition());
            }

            //- Synchronize values on selected edges. Edges are represented
            //  by the two vertices that make it up so global edges never get
            //  constructed.
            template<class T, class CombineOp>
            static void syncEdgeMap
            (
                const polyMesh& mesh,
                EdgeMap<T>& l,
                const CombineOp& cop
            )
            {
                syncEdgeMap(mesh, l, cop, mapDistribute::transform());
            }

            //- Synchronize locations on selected edges.
            template<class CombineOp>
            static void syncEdgePositions
            (
                const polyMesh& mesh,
                EdgeMap<point>& l,
                const CombineOp& cop
            )
            {
                syncEdgeMap(mesh, l, cop, mapDistribute::transformPosition());
            }

        // PackedList versions

            template<unsigned nBits, class CombineOp>
            static void syncFaceList
            (
                const polyMesh& mesh,
                PackedList<nBits>& faceValues,
                const CombineOp& cop
            );

            template<unsigned nBits>
            static void swapFaceList
            (
                const polyMesh& mesh,
                PackedList<nBits>& faceValues
            );

            template<unsigned nBits, class CombineOp>
            static void syncPointList
            (
                const polyMesh& mesh,
                PackedList<nBits>& pointValues,
                const CombineOp& cop,
                const unsigned int nullValue
            );

            template<unsigned nBits, class CombineOp>
            static void syncEdgeList
            (
                const polyMesh& mesh,
                PackedList<nBits>& edgeValues,
                const CombineOp& cop,
                const unsigned int nullValue
            );

        // Other

            //- Get per point whether it is uncoupled or a master of a
            //  coupled set of points
            static PackedBoolList getMasterPoints(const polyMesh&);

            //- Get per edge whether it is uncoupled or a master of a
            //  coupled set of edges
            static PackedBoolList getMasterEdges(const polyMesh&);

            //- Get per face whether it is uncoupled or a master of a
            //  coupled set of faces
            static PackedBoolList getMasterFaces(const polyMesh&);

            //- Get per face whether it is internal or a master of a
            //  coupled set of faces
            static PackedBoolList getInternalOrMasterFaces(const polyMesh&);

            //- Get per face whether it is internal or coupled
            static PackedBoolList getInternalOrCoupledFaces(const polyMesh&);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "syncToolsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
