/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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

Class
    Foam::isoSurface

Description
    A surface formed by the iso value.
    After "Regularised Marching Tetrahedra: improved iso-surface extraction",
    G.M. Treece, R.W. Prager and A.H. Gee.

    Note:
    - does tets without using cell centres/cell values. Not tested.
    - regularisation can create duplicate triangles/non-manifold surfaces.
    Current handling of those is bit ad-hoc for now and not perfect.
    - regularisation does not do boundary points so as to preserve the
      boundary perfectly.
    - uses geometric merge with fraction of bounding box as distance.
    - triangles can be between two cell centres so constant sampling
      does not make sense.
    - on empty patches behaves like zero gradient.
    - does not do 2D correctly, creates non-flat iso surface.
    - on processor boundaries might two overlapping (identical) triangles
      (one from either side)

    The handling on coupled patches is a bit complex. All fields
    (values and coordinates) get rewritten so
    - empty patches get zerogradient (value) and facecentre (coordinates)
    - separated processor patch faces get interpolate (value) and facecentre
      (coordinates). (this is already the default for cyclics)
    - non-separated processor patch faces keep opposite value and cell centre

    Now the triangle generation on non-separated processor patch faces
    can use the neighbouring value. Any separated processor face or cyclic
    face gets handled just like any boundary face.


SourceFiles
    isoSurface.C

\*---------------------------------------------------------------------------*/

#ifndef isoSurface_H
#define isoSurface_H

#include "triSurface.H"
#include "labelPair.H"
#include "pointIndexHit.H"
#include "PackedBoolList.H"
#include "volFields.H"
#include "slicedVolFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                       Class isoSurface Declaration
\*---------------------------------------------------------------------------*/

class isoSurface
:
    public triSurface
{
    // Private data

        enum segmentCutType
        {
            NEARFIRST,      // intersection close to e.first()
            NEARSECOND,     //  ,,                   e.second()
            NOTNEAR         // no intersection
        };

        enum cellCutType
        {
            NOTCUT,     // not cut
            SPHERE,     // all edges to cell centre cut
            CUT         // normal cut
        };


        //- Reference to mesh
        const fvMesh& mesh_;

        const scalarField& pVals_;

        //- Input volScalarField with separated coupled patches rewritten
        autoPtr<slicedVolScalarField> cValsPtr_;

        //- Isosurface value
        const scalar iso_;

        //- Regularise?
        const Switch regularise_;

        //- When to merge points
        const scalar mergeDistance_;


        //- Whether face might be cut
        List<cellCutType> faceCutType_;

        //- Whether cell might be cut
        List<cellCutType> cellCutType_;

        //- Estimated number of cut cells
        label nCutCells_;

        //- For every triangle the original cell in mesh
        labelList meshCells_;

        //- For every unmerged triangle point the point in the triSurface
        labelList triPointMergeMap_;


    // Private Member Functions

        // Point synchronisation

            //- Does tensor differ (to within mergeTolerance) from identity
            bool noTransform(const tensor& tt) const;

            //- Is patch a collocated (i.e. non-separated) coupled patch?
            static bool collocatedPatch(const polyPatch&);

            //- Per face whether is collocated
            PackedBoolList collocatedFaces(const coupledPolyPatch&) const;

            //- Synchronise points on all non-separated coupled patches
            void syncUnseparatedPoints
            (
                pointField& collapsedPoint,
                const point& nullValue
            ) const;


        //- Get location of iso value as fraction in between s0,s1
        scalar isoFraction
        (
            const scalar s0,
            const scalar s1
        ) const;

        //- Check if any edge of a face is cut
        bool isEdgeOfFaceCut
        (
            const scalarField& pVals,
            const face& f,
            const bool ownLower,
            const bool neiLower
        ) const;

        //- Get neighbour value and position.
        void getNeighbour
        (
            const labelList& boundaryRegion,
            const volVectorField& meshC,
            const volScalarField& cVals,
            const label celli,
            const label facei,
            scalar& nbrValue,
            point& nbrPoint
        ) const;

        //- Determine for every face/cell whether it (possibly) generates
        //  triangles.
        void calcCutTypes
        (
            const labelList& boundaryRegion,
            const volVectorField& meshC,
            const volScalarField& cVals,
            const scalarField& pVals
        );

        static point calcCentre(const triSurface&);

        //- Determine per cc whether all near cuts can be snapped to single
        //  point.
        void calcSnappedCc
        (
            const labelList& boundaryRegion,
            const volVectorField& meshC,
            const volScalarField& cVals,
            const scalarField& pVals,
            DynamicList<point>& snappedPoints,
            labelList& snappedCc
        ) const;

        //- Determine per point whether all near cuts can be snapped to single
        //  point.
        void calcSnappedPoint
        (
            const PackedBoolList& isBoundaryPoint,
            const labelList& boundaryRegion,
            const volVectorField& meshC,
            const volScalarField& cVals,
            const scalarField& pVals,
            DynamicList<point>& snappedPoints,
            labelList& snappedPoint
        ) const;


        //- Return input field with coupled (and empty) patch values rewritten
        template<class Type>
        tmp<SlicedGeometricField
        <Type, fvPatchField, slicedFvPatchField, volMesh>>
        adaptPatchFields
        (
            const GeometricField<Type, fvPatchField, volMesh>& fld
        ) const;

        //- Generate single point by interpolation or snapping
        template<class Type>
        Type generatePoint
        (
            const scalar s0,
            const Type& p0,
            const bool hasSnap0,
            const Type& snapP0,

            const scalar s1,
            const Type& p1,
            const bool hasSnap1,
            const Type& snapP1
        ) const;


        //- Note: cannot use simpler isoSurfaceCell::generateTriPoints since
        //  the need here to sometimes pass in remote 'snappedPoints'
        template<class Type>
        void generateTriPoints
        (
            const scalar s0,
            const Type& p0,
            const bool hasSnap0,
            const Type& snapP0,

            const scalar s1,
            const Type& p1,
            const bool hasSnap1,
            const Type& snapP1,

            const scalar s2,
            const Type& p2,
            const bool hasSnap2,
            const Type& snapP2,

            const scalar s3,
            const Type& p3,
            const bool hasSnap3,
            const Type& snapP3,

            DynamicList<Type>& points
        ) const;

        template<class Type>
        label generateFaceTriPoints
        (
            const volScalarField& cVals,
            const scalarField& pVals,

            const GeometricField<Type, fvPatchField, volMesh>& cCoords,
            const Field<Type>& pCoords,

            const DynamicList<Type>& snappedPoints,
            const labelList& snappedCc,
            const labelList& snappedPoint,
            const label facei,

            const scalar neiVal,
            const Type& neiPt,
            const bool hasNeiSnap,
            const Type& neiSnapPt,

            DynamicList<Type>& triPoints,
            DynamicList<label>& triMeshCells
        ) const;

        template<class Type>
        void generateTriPoints
        (
            const volScalarField& cVals,
            const scalarField& pVals,

            const GeometricField<Type, fvPatchField, volMesh>& cCoords,
            const Field<Type>& pCoords,

            const DynamicList<Type>& snappedPoints,
            const labelList& snappedCc,
            const labelList& snappedPoint,

            DynamicList<Type>& triPoints,
            DynamicList<label>& triMeshCells
        ) const;

        template<class Type>
        static tmp<Field<Type>> interpolate
        (
            const label nPoints,
            const labelList& triPointMergeMap,
            const DynamicList<Type>& unmergedValues
        );

        triSurface stitchTriPoints
        (
            const bool checkDuplicates,
            const List<point>& triPoints,
            labelList& triPointReverseMap,  // unmerged to merged point
            labelList& triMap               // merged to unmerged triangle
        ) const;

        //- Check single triangle for (topological) validity
        static bool validTri(const triSurface&, const label);

        static triSurface subsetMesh
        (
            const triSurface& s,
            const labelList& newToOldFaces,
            labelList& oldToNewPoints,
            labelList& newToOldPoints
        );

public:

    //- Declare friendship with isoSurfaceCell to share some functionality
    friend class isoSurfaceCell;


    //- Runtime type information
    TypeName("isoSurface");


    // Constructors

        //- Construct from cell values and point values. Uses boundaryField
        //  for boundary values. Holds reference to cellIsoVals and
        //  pointIsoVals.
        isoSurface
        (
            const volScalarField& cellIsoVals,
            const scalarField& pointIsoVals,
            const scalar iso,
            const bool regularise,
            const scalar mergeTol = 1e-6    // fraction of bounding box
        );


    // Member Functions

        //- For every triangle the original cell in mesh
        const labelList& meshCells() const
        {
            return meshCells_;
        }

        //- Interpolates cCoords,pCoords. Uses the references to the original
        //  fields used to create the iso surface.
        template<class Type>
        tmp<Field<Type>> interpolate
        (
            const GeometricField<Type, fvPatchField, volMesh>& cCoords,
            const Field<Type>& pCoords
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "isoSurfaceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
