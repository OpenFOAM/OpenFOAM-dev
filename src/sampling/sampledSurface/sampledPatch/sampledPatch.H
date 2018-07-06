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
    Foam::sampledPatch

Description
    A sampledSurface on patches. Non-triangulated by default.

SourceFiles
    sampledPatch.C

\*---------------------------------------------------------------------------*/

#ifndef sampledPatch_H
#define sampledPatch_H

#include "sampledSurface.H"
#include "MeshedSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class sampledPatch Declaration
\*---------------------------------------------------------------------------*/

class sampledPatch
:
    public MeshedSurface<face>,
    public sampledSurface
{
    //- Private typedefs for convenience
    typedef MeshedSurface<face> MeshStorage;


    // Private data

        //- Name of patches
        const wordReList patchNames_;

        //- Corresponding patchIDs
        mutable labelList patchIDs_;

        //- Triangulated faces or keep faces as is
        bool triangulate_;

        //- Track if the surface needs an update
        mutable bool needsUpdate_;

        //- For every face (or triangle) the originating patch
        labelList patchIndex_;

        //- For every face (or triangle) the index in the originating patch
        labelList patchFaceLabels_;

        //- Start indices (in patchFaceLabels_) of patches
        labelList patchStart_;


    // Private Member Functions

        //- Sample field on faces
        template<class Type>
        tmp<Field<Type>> sampleField
        (
            const GeometricField<Type, fvPatchField, volMesh>& vField
        ) const;

        //- Sample surface field on faces
        template<class Type>
        tmp<Field<Type>> sampleField
        (
            const GeometricField<Type, fvsPatchField, surfaceMesh>& sField
        ) const;

        template<class Type>
        tmp<Field<Type>> interpolateField(const interpolation<Type>&) const;

        //- Re-map action on triangulation or cleanup
        virtual void remapFaces(const labelUList& faceMap);


protected:

        const wordReList& patchNames() const
        {
            return patchNames_;
        }

        const labelList& patchIDs() const;

        const labelList& patchStart() const
        {
            return patchStart_;
        }

        const labelList& patchFaceLabels() const
        {
            return patchFaceLabels_;
        }


public:

    //- Runtime type information
    TypeName("sampledPatch");


    // Constructors

        //- Construct from components
        sampledPatch
        (
            const word& name,
            const polyMesh& mesh,
            const wordReList& patchNames,
            const bool triangulate = false
        );

        //- Construct from dictionary
        sampledPatch
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~sampledPatch();


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
            return MeshStorage::points();
        }

        //- Faces of surface
        virtual const faceList& faces() const
        {
            return MeshStorage::faces();
        }


        // Sample
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

            //- Surface sample field on surface
            virtual tmp<scalarField> sample
            (
                const surfaceScalarField&
            ) const;

            //- Surface Sample field on surface
            virtual tmp<vectorField> sample
            (
                const surfaceVectorField&
            ) const;

            //- Surface sample field on surface
            virtual tmp<sphericalTensorField> sample
            (
                const surfaceSphericalTensorField&
            ) const;

            //- Surface sample field on surface
            virtual tmp<symmTensorField> sample
            (
                const surfaceSymmTensorField&
            ) const;

            //- Surface sample field on surface
            virtual tmp<tensorField> sample
            (
                const surfaceTensorField&
            ) const;


        // Interpolate

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

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "sampledPatchTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
