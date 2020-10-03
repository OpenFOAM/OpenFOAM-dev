/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::surfaceSlipDisplacementPointPatchVectorField

Description
    Displacement follows a triSurface. Use in a displacementMotionSolver
    as a bc on the pointDisplacement field.
    Following is done by calculating the projection onto the surface according
    to the projectMode
    - NEAREST : nearest
    - POINTNORMAL : intersection with point normal
    - FIXEDNORMAL : intersection with fixed vector

    Optionally (intersection only) removes a component ("wedgePlane") to
    stay in 2D.

    Needs:
    - geometry : dictionary with searchableSurfaces. (usually
      triSurfaceMeshes in constant/triSurface)
    - projectMode : see above
    - projectDirection : if projectMode = fixedNormal
    - wedgePlane : -1 or component to knock out of intersection normal
    - frozenPointsZone : empty or name of pointZone containing points
                         that do not move

SourceFiles
    surfaceSlipDisplacementPointPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceSlipDisplacementPointPatchVectorField_H
#define surfaceSlipDisplacementPointPatchVectorField_H

#include "pointPatchFields.H"
#include "searchableSurfaces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class surfaceSlipDisplacementPointPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class surfaceSlipDisplacementPointPatchVectorField
:
    public pointPatchVectorField
{
public:

    // Public data types

        enum projectMode
        {
            NEAREST,
            POINTNORMAL,
            FIXEDNORMAL
        };

private:

    // Private Data

        //- Project mode names
        static const NamedEnum<projectMode, 3> projectModeNames_;

        //- Names of surfaces
        const dictionary surfacesDict_;

        //- How to project/project onto surface
        const projectMode projectMode_;

        //- Direction to project
        const vector projectDir_;

        //- Plane for 2D wedge case or -1.
        const label wedgePlane_;

        //- pointZone with frozen points
        const word frozenPointsZone_;

        //- Demand driven: surface to project
        mutable autoPtr<searchableSurfaces> surfacesPtr_;


    // Private Member Functions

        //- Calculate displacement (w.r.t. points0()) to project onto surface
        void calcProjection(vectorField& displacement) const;


public:

    //- Runtime type information
    TypeName("surfaceSlipDisplacement");


    // Constructors

        //- Construct from patch and internal field
        surfaceSlipDisplacementPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        surfaceSlipDisplacementPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<vector> onto a new patch
        surfaceSlipDisplacementPointPatchVectorField
        (
            const surfaceSlipDisplacementPointPatchVectorField&,
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        surfaceSlipDisplacementPointPatchVectorField
        (
            const surfaceSlipDisplacementPointPatchVectorField&
        ) = delete;

        //- Copy constructor setting internal field reference
        surfaceSlipDisplacementPointPatchVectorField
        (
            const surfaceSlipDisplacementPointPatchVectorField&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<pointPatchVectorField> clone
        (
            const DimensionedField<vector, pointMesh>& iF
        ) const
        {
            return autoPtr<pointPatchVectorField>
            (
                new surfaceSlipDisplacementPointPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }

    // Member Functions

        //- Surface to follow. Demand loads surfaceNames.
        const searchableSurfaces& surfaces() const;

        //- Update the patch field
        virtual void evaluate
        (
            const Pstream::commsTypes commsType=Pstream::commsTypes::blocking
        );

        //- Write
        virtual void write(Ostream&) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=
        (
            const surfaceSlipDisplacementPointPatchVectorField&
        ) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
