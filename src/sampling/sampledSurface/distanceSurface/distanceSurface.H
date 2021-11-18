/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
     Foam::sampledSurfaces::distanceSurface

Description
    A sampledSurface defined by a distance to a surface.

    Example:
    \verbatim
    {
        type        distanceSurface;
        surfaceType searchableBox;
        min         (-1 -1 -1);
        max         (1 1 1);
        distance    0.1;
        signed      yes;
        filtering   full;
        interpolate yes;
    }
    \endverbatim

Usage
    \table
        Property     | Description             | Required    | Default value
        surfaceType  | the type of surface to sample from | yes |
        distance     | the distance from which to sample the surface | yes |
        signed       | sample only on one side of the surface, specified by \\
                       the sign of the distance | yes |
        filtering    | the level of filtering to perform on the iso-surface \\
                                               | no          | full
        interpolate  | interpolate values to the surface points | no | no
    \endtable

See also
    Foam::isoSurface

SourceFiles
    distanceSurface.C

\*---------------------------------------------------------------------------*/

#ifndef sampledSurfaces_distanceSurface_H
#define sampledSurfaces_distanceSurface_H

#include "sampledSurface.H"
#include "searchableSurface.H"
#include "isoSurface.H"
#include "fvMeshSubset.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace sampledSurfaces
{

/*---------------------------------------------------------------------------*\
                       Class distanceSurface Declaration
\*---------------------------------------------------------------------------*/

class distanceSurface
:
    public sampledSurface
{
    // Private Data

        //- Surface
        const autoPtr<searchableSurface> surfPtr_;

        //- Distance value
        const scalar distance_;

        //- Signed distance
        const bool signed_;

        //- Whether to coarsen
        const isoSurface::filterType filter_;

        //- Whether to recalculate cell values as average of point values
        const Switch average_;

        //- If restricted to zones, name of this zone or a regular expression
        wordRe zoneKey_;

        //- For zones: patch to put exposed faces into
        mutable word exposedPatchName_;

        //- Track if the surface needs an update
        mutable bool needsUpdate_;

        //- Optional subsetted mesh
        autoPtr<fvMeshSubset> subMeshPtr_;

        //- Distance to cell centres
        autoPtr<volScalarField> cellDistancePtr_;

        //- Distance to points
        scalarField pointDistance_;

        //- Constructed iso surface
        autoPtr<isoSurface> isoSurfPtr_;


    // Private Member Functions

        //- Create iso surface
        void createGeometry();

        //- Sample field on faces
        template<class Type>
        tmp<Field<Type>> sampleField
        (
            const GeometricField<Type, fvPatchField, volMesh>& vField
        ) const;


        template<class Type>
        tmp<Field<Type>>
        interpolateField(const interpolation<Type>&) const;


public:

    //- Runtime type information
    TypeName("distanceSurface");


    // Constructors

        //- Construct from dictionary
        distanceSurface
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~distanceSurface();


    // Member Functions

        //- Does the surface need an update?
        virtual bool needsUpdate() const;

        //- Mark the surface as needing an update.
        //  May also free up unneeded data.
        //  Return false if surface was already marked as expired.
        virtual bool expire();

        //- Update the surface as required.
        //  Do nothing (and return false) if no update was needed
        virtual bool update();

        //- Points of surface
        virtual const pointField& points() const
        {
            return surface().points();
        }

        //- Faces of surface
        virtual const faceList& faces() const
        {
            return surface().faces();
        }

        const isoSurface& surface() const
        {
            return isoSurfPtr_();
        }

        //- Sample field on surface
        virtual tmp<scalarField> sample
        (
            const volScalarField&
        ) const;

        //- Sample field on surface
        virtual tmp<vectorField> sample
        (
            const volVectorField&
        ) const;

        //- Sample field on surface
        virtual tmp<sphericalTensorField> sample
        (
            const volSphericalTensorField&
        ) const;

        //- Sample field on surface
        virtual tmp<symmTensorField> sample
        (
            const volSymmTensorField&
        ) const;

        //- Sample field on surface
        virtual tmp<tensorField> sample
        (
            const volTensorField&
        ) const;


        //- Interpolate field on surface
        virtual tmp<scalarField> interpolate
        (
            const interpolation<scalar>&
        ) const;

        //- Interpolate field on surface
        virtual tmp<vectorField> interpolate
        (
            const interpolation<vector>&
        ) const;

        //- Interpolate field on surface
        virtual tmp<sphericalTensorField> interpolate
        (
            const interpolation<sphericalTensor>&
        ) const;

        //- Interpolate field on surface
        virtual tmp<symmTensorField> interpolate
        (
            const interpolation<symmTensor>&
        ) const;

        //- Interpolate field on surface
        virtual tmp<tensorField> interpolate
        (
            const interpolation<tensor>&
        ) const;

        //- Write
        virtual void print(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace sampledSurfaces
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "distanceSurfaceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
