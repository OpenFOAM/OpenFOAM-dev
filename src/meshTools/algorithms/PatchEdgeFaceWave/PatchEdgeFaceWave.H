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

Class
    Foam::PatchEdgeFaceWave

Description
    Wave propagation of information along patch. Every iteration
    information goes through one layer of faces. Templated on information
    that is transferred.

SourceFiles
    PatchEdgeFaceWave.C

\*---------------------------------------------------------------------------*/

#ifndef PatchEdgeFaceWave_H
#define PatchEdgeFaceWave_H

#include "scalarField.H"
#include "PackedBoolList.H"
#include "PrimitivePatch.H"
#include "vectorTensorTransform.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;

/*---------------------------------------------------------------------------*\
                    Class PatchEdgeFaceWaveName Declaration
\*---------------------------------------------------------------------------*/

TemplateName(PatchEdgeFaceWave);


/*---------------------------------------------------------------------------*\
                           Class PatchEdgeFaceWave Declaration
\*---------------------------------------------------------------------------*/

template
<
    class PrimitivePatchType,
    class Type,
    class TrackingData = int
>
class PatchEdgeFaceWave
:
    public PatchEdgeFaceWaveName
{
  // Private static data

        //- Relative tolerance. Stop propagation if relative changes
        //  less than this tolerance (responsability for checking this is
        //  up to Type implementation)
        static scalar propagationTol_;

        //- Used as default trackdata value to satisfy default template
        //  argument.
        static label dummyTrackData_;


    // Private data

        //- Reference to mesh
        const polyMesh& mesh_;

        //- Reference to patch
        const PrimitivePatchType& patch_;

        //- Wall information for all edges
        UList<Type>& allEdgeInfo_;

        //- Information on all patch faces
        UList<Type>& allFaceInfo_;

        //- Additional data to be passed into container
        TrackingData& td_;

        //- Has edge changed
        PackedBoolList changedEdge_;

        //- List of changed edges
        DynamicList<label> changedEdges_;

        //- Has face changed
        PackedBoolList changedFace_;

        //- List of changed faces
        DynamicList<label> changedFaces_;

        //- Number of evaluations
        label nEvals_;

        //- Number of unvisited faces/edges
        label nUnvisitedEdges_;
        label nUnvisitedFaces_;


        // Addressing between edges of patch_ and globalData.coupledPatch()
        labelList patchEdges_;
        labelList coupledEdges_;
        PackedBoolList sameEdgeOrientation_;


    // Private Member Functions

        //- Updates edgeInfo with information from neighbour. Updates all
        //  statistics.
        bool updateEdge
        (
            const label edgeI,
            const label neighbourFaceI,
            const Type& neighbourInfo,
            Type& edgeInfo
        );

        //- Updates faceInfo with information from neighbour. Updates all
        //  statistics.
        bool updateFace
        (
            const label faceI,
            const label neighbourEdgeI,
            const Type& neighbourInfo,
            Type& faceInfo
        );

        //- Update coupled edges
        void syncEdges();

        //- Disallow default bitwise copy construct
        PatchEdgeFaceWave(const PatchEdgeFaceWave&);

        //- Disallow default bitwise assignment
        void operator=(const PatchEdgeFaceWave&);


public:

    // Static Functions

        //- Access to tolerance
        static scalar propagationTol()
        {
            return propagationTol_;
        }

        //- Change tolerance
        static void setPropagationTol(const scalar tol)
        {
            propagationTol_ = tol;
        }


    // Constructors

        //- Construct from patch, list of changed edges with the Type
        //  for these edges. Gets work arrays to operate on, one of size
        //  number of patch edges, the other number of patch faces.
        //  Iterates until nothing changes or maxIter reached.
        //  (maxIter can be 0)
        PatchEdgeFaceWave
        (
            const polyMesh& mesh,
            const PrimitivePatchType& patch,
            const labelList& initialEdges,
            const List<Type>& initialEdgesInfo,
            UList<Type>& allEdgeInfo,
            UList<Type>& allFaceInfo,
            const label maxIter,
            TrackingData& td = dummyTrackData_
        );

        //- Construct from patch. Use setEdgeInfo and iterate() to do
        //  actual calculation
        PatchEdgeFaceWave
        (
            const polyMesh& mesh,
            const PrimitivePatchType& patch,
            UList<Type>& allEdgeInfo,
            UList<Type>& allFaceInfo,
            TrackingData& td = dummyTrackData_
        );


    // Member Functions

        //- Access allEdgeInfo
        UList<Type>& allEdgeInfo() const
        {
            return allEdgeInfo_;
        }

        //- Access allFaceInfo
        UList<Type>& allFaceInfo() const
        {
            return allFaceInfo_;
        }

        //- Additional data to be passed into container
        const TrackingData& data() const
        {
            return td_;
        }

        //- Get number of unvisited faces, i.e. faces that were not (yet)
        //  reached from walking across patch. This can happen from
        //  - not enough iterations done
        //  - a disconnected patch
        //  - a patch without walls in it
        label getUnsetFaces() const;

        label getUnsetEdges() const;

        //- Copy initial data into allEdgeInfo_
        void setEdgeInfo
        (
            const labelList& changedEdges,
            const List<Type>& changedEdgesInfo
        );

        //- Propagate from edge to face. Returns total number of faces
        //  (over all processors) changed.
        label edgeToFace();

        //- Propagate from face to edge. Returns total number of edges
        //  (over all processors) changed.
        label faceToEdge();

        //- Iterate until no changes or maxIter reached. Returns actual
        //  number of iterations.
        label iterate(const label maxIter);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Update operation
template
<
    class PrimitivePatchType,
    class Type,
    class TrackingData = int
>
class updateOp
{
    //- Additional data to be passed into container
    const polyMesh& mesh_;
    const PrimitivePatchType& patch_;
    const scalar tol_;
    TrackingData& td_;

public:
    updateOp
    (
        const polyMesh& mesh,
        const PrimitivePatchType& patch,
        const scalar tol,
        TrackingData& td
    )
    :
        mesh_(mesh),
        patch_(patch),
        tol_(tol),
        td_(td)
    {}

    void operator()(Type& x, const Type& y) const
    {
        if (y.valid(td_))
        {
            x.updateEdge(mesh_, patch_, y, true, tol_, td_);
        }
    }
};


//- Transform operation
template
<
    class PrimitivePatchType,
    class Type,
    class TrackingData = int
>
class transformOp
{
    //- Additional data to be passed into container
    const polyMesh& mesh_;
    const PrimitivePatchType& patch_;
    const scalar tol_;
    TrackingData& td_;

public:
    transformOp
    (
        const polyMesh& mesh,
        const PrimitivePatchType& patch,
        const scalar tol,
        TrackingData& td
    )
    :
        mesh_(mesh),
        patch_(patch),
        tol_(tol),
        td_(td)
    {}

    void operator()
    (
        const vectorTensorTransform& vt,
        const bool forward,
        List<Type>& fld
    ) const
    {
        if (forward)
        {
            forAll(fld, i)
            {
                fld[i].transform(mesh_, patch_, vt.R(), tol_, td_);
            }
        }
        else
        {
            forAll(fld, i)
            {
                fld[i].transform(mesh_, patch_, vt.R().T(), tol_, td_);
            }
        }
    }
};

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "PatchEdgeFaceWave.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
