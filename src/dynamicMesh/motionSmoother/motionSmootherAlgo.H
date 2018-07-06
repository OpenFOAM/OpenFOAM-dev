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

Class
    Foam::motionSmootherAlgo

Description
    Given a displacement moves the mesh by scaling the displacement back
    until there are no more mesh errors.

    Holds displacement field (read upon construction since need boundary
    conditions) and scaling factor and optional patch number on which to
    scale back displacement.

    E.g.
    \verbatim
        // Construct iterative mesh mover.
        motionSmoother meshMover(mesh, labelList(1, patchi));

        // Set desired displacement:
        meshMover.displacement() = ..

        for (label iter = 0; iter < maxIter; iter++)
        {
            if (meshMover.scaleMesh(true))
            {
                Info<< "Successfully moved mesh" << endl;
                return true;
            }
        }
    \endverbatim

Note
    - Shared points (parallel): a processor can have points which are part of
    pp on another processor but have no pp itself (i.e. it has points
    and/or edges but no faces of pp). Hence we have to be careful when e.g.
    synchronising displacements that the value from the processor which has
    faces of pp get priority. This is currently handled in setDisplacement
    by resetting the internal displacement to zero before doing anything
    else. The combine operator used will give preference to non-zero
    values.

    - Various routines take baffles. These are sets of boundary faces that
    are treated as a single internal face. This is a hack used to apply
    movement to internal faces.

    - Mesh constraints are looked up from the supplied dictionary. (uses
    recursive lookup)

SourceFiles
    motionSmootherAlgo.C
    motionSmootherAlgoTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef motionSmootherAlgo_H
#define motionSmootherAlgo_H

#include "pointFields.H"
#include "HashSet.H"
#include "PackedBoolList.H"
#include "indirectPrimitivePatch.H"
#include "className.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMeshGeometry;
class faceSet;

/*---------------------------------------------------------------------------*\
                           Class motionSmootherAlgo Declaration
\*---------------------------------------------------------------------------*/

class motionSmootherAlgo
{
    // Private class

        //- To synchronise displacements. We want max displacement since
        //  this is what is specified on pp and internal mesh will have
        //  zero displacement.
        class maxMagEqOp
        {

        public:

            void operator()(vector& x, const vector& y) const
            {
                for (direction i = 0; i < vector::nComponents; i++)
                {
                    scalar magX = mag(x[i]);
                    scalar magY = mag(y[i]);

                    if (magX < magY)
                    {
                        x[i] = y[i];
                    }
                    else if (magX == magY)
                    {
                        if (y[i] > x[i])
                        {
                            x[i] = y[i];
                        }
                    }
                }
            }
        };


    // Private data

        //- Reference to polyMesh. Non-const since we move mesh.
        polyMesh& mesh_;

        //- Reference to pointMesh
        pointMesh& pMesh_;

        //- Reference to face subset of all adaptPatchIDs
        indirectPrimitivePatch& pp_;

        //- Displacement field
        pointVectorField& displacement_;

        //- Scale factor for displacement
        pointScalarField& scale_;

        //- Starting mesh position
        pointField& oldPoints_;


        // Internal data

            //- Indices of fixedValue patches that we're allowed to modify the
            // displacement on.
            const labelList adaptPatchIDs_;

            // Smoothing and checking parameters
            dictionary paramDict_;

            //- Is mesh point on boundary or not
            PackedBoolList isInternalPoint_;

            //- Is edge master (always except if on coupled boundary and on
            //  lower processor)
            PackedBoolList isMasterEdge_;


    // Private Member Functions

        //- Average of connected points.
        template<class Type>
        tmp<GeometricField<Type, pointPatchField, pointMesh>> avg
        (
            const GeometricField<Type, pointPatchField, pointMesh>& fld,
            const scalarField& edgeWeight
        ) const;

        //- Average position of connected points.
        template<class Type>
        tmp<GeometricField<Type, pointPatchField, pointMesh>> avgPositions
        (
            const GeometricField<Type, pointPatchField, pointMesh>& fld,
            const scalarField& edgeWeight
        ) const;

        //- Check constraints
        template<class Type>
        static void checkConstraints
        (
            GeometricField<Type, pointPatchField, pointMesh>&
        );

        //- Test synchronisation of generic field (not positions!) on points
        template<class Type, class CombineOp>
        void testSyncField
        (
            const Field<Type>&,
            const CombineOp& cop,
            const Type& zero,
            const scalar maxMag
        ) const;

        //- Test synchronisation of points
        void testSyncPositions(const pointField&, const scalar maxMag) const;

        static void checkFld(const pointScalarField&);

        //- Get points used by given faces
        labelHashSet getPoints(const labelHashSet&) const;

        //- Calculate per-edge weight
        tmp<scalarField> calcEdgeWeights(const pointField&) const;

        //- Explicit smoothing and min on all affected internal points
        void minSmooth
        (
            const scalarField& edgeWeights,
            const PackedBoolList& isAffectedPoint,
            const pointScalarField& fld,
            pointScalarField& newFld
        ) const;

        //- Same but only on selected points (usually patch points)
        void minSmooth
        (
            const scalarField& edgeWeights,
            const PackedBoolList& isAffectedPoint,
            const labelList& meshPoints,
            const pointScalarField& fld,
            pointScalarField& newFld
        ) const;

        //- Scale certain (internal) points of a field
        void scaleField
        (
            const labelHashSet& pointLabels,
            const scalar scale,
            pointScalarField&
        ) const;

        //- As above but points have to be in meshPoints as well
        //  (usually to scale patch points)
        void scaleField
        (
            const labelList& meshPoints,
            const labelHashSet& pointLabels,
            const scalar scale,
            pointScalarField&
        ) const;

        //- Lower on internal points
        void subtractField
        (
            const labelHashSet& pointLabels,
            const scalar f,
            pointScalarField&
        ) const;

        //- As above but points have to be in meshPoints as well
        //  (usually to scale patch points)
        void subtractField
        (
            const labelList& meshPoints,
            const labelHashSet& pointLabels,
            const scalar scale,
            pointScalarField&
        ) const;

        //- Helper function. Is point internal?
        bool isInternalPoint(const label pointi) const;

        //- Given a set of faces that cause smoothing and a number of
        //  iterations determine the maximum set of points who are affected
        //  and the accordingly affected faces.
        void getAffectedFacesAndPoints
        (
            const label nPointIter,
            const faceSet& wrongFaces,

            labelList& affectedFaces,
            PackedBoolList& isAffectedPoint
        ) const;

        //- Disallow default bitwise copy construct
        motionSmootherAlgo(const motionSmootherAlgo&);

        //- Disallow default bitwise assignment
        void operator=(const motionSmootherAlgo&);


public:

    ClassName("motionSmootherAlgo");

    // Constructors

        //- Construct from mesh, patches to work on and smoothing parameters.
        motionSmootherAlgo
        (
            polyMesh&,
            pointMesh&,
            indirectPrimitivePatch& pp,         // 'outside' points
            pointVectorField& displacement,
            pointScalarField& scale,
            pointField& oldPoints,
            const labelList& adaptPatchIDs,     // patches forming 'outside'
            const dictionary& paramDict
        );


    //- Destructor
    ~motionSmootherAlgo();


    // Member Functions

        // Access

            //- Reference to mesh
            const polyMesh& mesh() const;

            //- Reference to pointMesh
            const pointMesh& pMesh() const;

            //- Reference to patch
            const indirectPrimitivePatch& patch() const;

            //- Patch labels that are being adapted
            const labelList& adaptPatchIDs() const;

            const dictionary& paramDict() const;



        // Edit

            //- Take over existing mesh position.
            void correct();


            //- Set patch fields on patchIDs to be consistent with
            //  all other boundary conditions
            static void setDisplacementPatchFields
            (
                const labelList& patchIDs,
                pointVectorField& pointDisplacement
            );

            //- Set patch fields on displacement to be consistent with
            //  internal values.
            void setDisplacementPatchFields();

            //- Set displacement field from displacement on patch points.
            //  Modify provided displacement to be consistent with actual
            //  boundary conditions on displacement. Note: resets the
            //  displacement to be 0 on coupled patches beforehand
            //  to make sure shared points
            //  partially on pp (on some processors) and partially not
            //  (on other processors) get the value from pp.
            static void setDisplacement
            (
                const labelList& patchIDs,
                const indirectPrimitivePatch& pp,
                pointField& patchDisp,
                pointVectorField& displacement
            );

            void setDisplacement(pointField& patchDisp);

            //- Special correctBoundaryConditions which evaluates fixedValue
            //  patches first so they get overwritten with any constraint
            //  bc's.
            void correctBoundaryConditions(pointVectorField&) const;

            //- Apply optional point constraint (2d correction)
            void modifyMotionPoints(pointField& newPoints) const;

            //- Get the current points (oldPoints+scale*displacement)
            tmp<pointField> curPoints() const;

            //- Set the errorReduction (by how much to scale the displacement
            //  at error locations) parameter. Returns the old value.
            //  Set to 0 (so revert to old mesh) grows out one cell layer
            //  from error faces.
            scalar setErrorReduction(const scalar);

            //- Move mesh with given scale. Return true if mesh ok or has
            //  less than nAllow errors, false
            //  otherwise and locally update scale. Smoothmesh=false means only
            //  patch points get moved.
            //  Parallel ok (as long as displacement field is consistent
            //  across patches)
            bool scaleMesh
            (
                labelList& checkFaces,
                const bool smoothMesh = true,
                const label nAllow = 0
            );

            //- Move mesh (with baffles) with given scale.
            bool scaleMesh
            (
                labelList& checkFaces,
                const List<labelPair>& baffles,
                const bool smoothMesh = true,
                const label nAllow = 0
            );

            //- Move mesh with externally provided mesh constraints
            bool scaleMesh
            (
                labelList& checkFaces,
                const List<labelPair>& baffles,
                const dictionary& paramDict,
                const dictionary& meshQualityDict,
                const bool smoothMesh = true,
                const label nAllow = 0
            );


            //- Update for new mesh geometry
            void movePoints();

            //- Update for new mesh topology
            void updateMesh();


            //- Check mesh with mesh settings in dict. Collects incorrect faces
            //  in set. Returns true if one or more faces in error.
            //  Parallel ok.
            static bool checkMesh
            (
                const bool report,
                const polyMesh& mesh,
                const dictionary& dict,
                labelHashSet& wrongFaces
            );

            //- Check (subset of mesh) with mesh settings in dict.
            //  Collects incorrect faces in set. Returns true if one
            //  or more faces in error. Parallel ok.
            static bool checkMesh
            (
                const bool report,
                const polyMesh& mesh,
                const dictionary& dict,
                const labelList& checkFaces,
                labelHashSet& wrongFaces
            );

            //- Check (subset of mesh including baffles) with mesh settings
            //  in dict. Collects incorrect faces in set. Returns true if one
            //  or more faces in error. Parallel ok.
            static bool checkMesh
            (
                const bool report,
                const polyMesh& mesh,
                const dictionary& dict,
                const labelList& checkFaces,
                const List<labelPair>& baffles,
                labelHashSet& wrongFaces
            );

            //- Check part of mesh with mesh settings in dict.
            //  Collects incorrect faces in set. Returns true if one or
            //  more faces in error. Parallel ok.
            static bool checkMesh
            (
                const bool report,
                const dictionary& dict,
                const polyMeshGeometry&,
                const labelList& checkFaces,
                labelHashSet& wrongFaces
            );

            //- Check part of mesh including baffles with mesh settings in dict.
            //  Collects incorrect faces in set. Returns true if one or
            //  more faces in error. Parallel ok.
            static bool checkMesh
            (
                const bool report,
                const dictionary& dict,
                const polyMeshGeometry&,
                const labelList& checkFaces,
                const List<labelPair>& baffles,
                labelHashSet& wrongFaces
            );

            // Helper functions to manipulate displacement vector.

                //- Fully explicit smoothing of fields (not positions)
                //  of internal points with varying diffusivity.
                template<class Type>
                void smooth
                (
                    const GeometricField<Type, pointPatchField, pointMesh>& fld,
                    const scalarField& edgeWeight,
                    GeometricField<Type, pointPatchField, pointMesh>& newFld
                ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "motionSmootherAlgoTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
