/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2022 OpenFOAM Foundation
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
    Foam::pointConstraints

Description
    Application of (multi-)patch point constraints.

    Note: includes all points which are on the boundary of a patch
          with a constraint. It includes them (even though the constraint
          will already be implemented through the patch evaluation)
          since these points might be
          coupled to points which are not on any constraint patch and we
          don't want to get inconsistency between the two points.

SourceFiles
    pointConstraints.C
    pointConstraintsTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef pointConstraints_H
#define pointConstraints_H

#include "DemandDrivenMeshObject.H"
#include "tensorField.H"
#include "pointFieldsFwd.H"
#include "pointConstraint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;

/*---------------------------------------------------------------------------*\
                      Class pointConstraints Declaration
\*---------------------------------------------------------------------------*/

class pointConstraints
:
    public DemandDrivenMeshObject
    <
        pointMesh,
        UpdateableMeshObject,
        pointConstraints
    >
{
    // Private Data

        //- Mesh points on which to apply special constraints
        labelList patchPatchPointConstraintPoints_;

        //- Special constraints (as tensors)
        tensorField patchPatchPointConstraintTensors_;

        //- Special constraints (raw)
        List<pointConstraint> patchPatchPointConstraints_;


    // Private Member Functions

        //- Make patch-patch constraints
        void makePatchPatchAddressing();


protected:

    friend class DemandDrivenMeshObject
    <
        pointMesh,
        UpdateableMeshObject,
        pointConstraints
    >;

    // Protected Constructors

        //- Constructor from pointMesh.
        explicit pointConstraints(const pointMesh&);


public:

    // Declare name of the class and its debug switch
    ClassName("pointConstraints");


    // Constructors

        //- Disallow default bitwise copy construction
        pointConstraints(const pointConstraints&) = delete;


    //- Destructor
    ~pointConstraints();


    // Member Functions

        // Access

            //- Mesh points on which to apply special constraints
            const labelList& patchPatchPointConstraintPoints() const
            {
                return patchPatchPointConstraintPoints_;
            }

            //- Special constraints
            const tensorField& patchPatchPointConstraintTensors() const
            {
                return patchPatchPointConstraintTensors_;
            }

            //- Actual constraints
            const List<pointConstraint>& patchPatchPointConstraints() const
            {
                return patchPatchPointConstraints_;
            }


        // Edit

            //- Correct weighting factors for moving mesh.
            virtual bool movePoints();

            //- Update mesh topology using the morph engine
            virtual void topoChange(const polyTopoChangeMap&);

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap&);

            //- Update corresponding to the given distribution map
            virtual void distribute(const polyDistributionMap&);


        // Interpolation functions

            //- Helper: sync data on collocated points only
            template<class Type, class CombineOp>
            static void syncUntransformedData
            (
                const polyMesh& mesh,
                List<Type>& pointData,
                const CombineOp& cop
            );

            //- Helper: set patchField values from internal values (on
            //  valuePointPatchFields). Opposite of
            //  pointPatchField::setInternalField
            template<class Type>
            static void setPatchFields
            (
                PointField<Type>&
            );

            //- Apply patch-patch constraints only
            template<class Type>
            void constrainCorners
            (
                PointField<Type>& pf
            ) const;

            //- Apply boundary conditions (single-patch constraints) and
            //  patch-patch constraints on generic pointField.
            //  Optionally override valuePointPatchFields with constrained
            //  values (see setPatchFields above)
            template<class Type>
            void constrain
            (
                PointField<Type>& pf,
                const bool overrideValue = false
            ) const;

            //- Apply boundary conditions (single-patch constraints),
            //  patch-patch constraints and
            //  two-D constraints on displacement field
            void constrainDisplacement
            (
                pointVectorField& displacement,
                const bool overrideValue = false
            ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const pointConstraints&) = delete;
};


template<>
void pointConstraints::constrainCorners<scalar>
(
    PointField<scalar>& pf
) const;

template<>
void pointConstraints::constrainCorners<label>
(
    PointField<label>& pf
) const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "pointConstraintsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
