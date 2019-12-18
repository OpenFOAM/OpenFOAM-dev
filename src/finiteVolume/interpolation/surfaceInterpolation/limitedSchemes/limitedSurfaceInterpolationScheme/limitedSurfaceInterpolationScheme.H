/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::limitedSurfaceInterpolationScheme

Description
    Abstract base class for limited surface interpolation schemes.

SourceFiles
    limitedSurfaceInterpolationScheme.C

\*---------------------------------------------------------------------------*/

#ifndef limitedSurfaceInterpolationScheme_H
#define limitedSurfaceInterpolationScheme_H

#include "surfaceInterpolationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class limitedSurfaceInterpolationScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class limitedSurfaceInterpolationScheme
:
    public surfaceInterpolationScheme<Type>
{

protected:

    // Protected data

        const surfaceScalarField& faceFlux_;


public:

    //- Runtime type information
    TypeName("limitedSurfaceInterpolationScheme");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            tmp,
            limitedSurfaceInterpolationScheme,
            Mesh,
            (
                const fvMesh& mesh,
                Istream& schemeData
            ),
            (mesh, schemeData)
        );

        declareRunTimeSelectionTable
        (
            tmp,
            limitedSurfaceInterpolationScheme,
            MeshFlux,
            (
                const fvMesh& mesh,
                const surfaceScalarField& faceFlux,
                Istream& schemeData
            ),
            (mesh, faceFlux, schemeData)
        );


    // Constructors

        //- Construct from mesh and faceFlux
        limitedSurfaceInterpolationScheme
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            faceFlux_(faceFlux)
        {}

        //- Construct from mesh and Istream.
        //  The name of the flux field is read from the Istream and looked-up
        //  from the mesh objectRegistry
        limitedSurfaceInterpolationScheme
        (
            const fvMesh& mesh,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            faceFlux_
            (
                mesh.lookupObject<surfaceScalarField>
                (
                    word(is)
                )
            )
        {}

        //- Disallow default bitwise copy construction
        limitedSurfaceInterpolationScheme
        (
            const limitedSurfaceInterpolationScheme&
        );


    // Selectors

        //- Return new tmp interpolation scheme
        static tmp<limitedSurfaceInterpolationScheme<Type>> New
        (
            const fvMesh& mesh,
            Istream& schemeData
        );

        //- Return new tmp interpolation scheme
        static tmp<limitedSurfaceInterpolationScheme<Type>> New
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& schemeData
        );


    //- Destructor
    virtual ~limitedSurfaceInterpolationScheme();


    // Member Functions

        //- Return the interpolation weighting factors
        virtual tmp<surfaceScalarField> limiter
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const = 0;

        //- Return the interpolation weighting factors for the given field,
        //  by limiting the given weights with the given limiter
        tmp<surfaceScalarField> weights
        (
            const GeometricField<Type, fvPatchField, volMesh>&,
            const surfaceScalarField& CDweights,
            tmp<surfaceScalarField> tLimiter
        ) const;

        //- Return the interpolation weighting factors for the given field
        virtual tmp<surfaceScalarField> weights
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        //- Return the interpolation weighting factors
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        flux
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const limitedSurfaceInterpolationScheme&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables

#define makelimitedSurfaceInterpolationTypeScheme(SS, Type)                    \
                                                                               \
defineNamedTemplateTypeNameAndDebug(SS<Type>, 0);                              \
                                                                               \
surfaceInterpolationScheme<Type>::addMeshConstructorToTable<SS<Type>>          \
    add##SS##Type##MeshConstructorToTable_;                                    \
                                                                               \
surfaceInterpolationScheme<Type>::addMeshFluxConstructorToTable<SS<Type>>      \
    add##SS##Type##MeshFluxConstructorToTable_;                                \
                                                                               \
limitedSurfaceInterpolationScheme<Type>::addMeshConstructorToTable<SS<Type>>   \
    add##SS##Type##MeshConstructorToLimitedTable_;                             \
                                                                               \
limitedSurfaceInterpolationScheme<Type>::                                      \
    addMeshFluxConstructorToTable<SS<Type>>                                    \
    add##SS##Type##MeshFluxConstructorToLimitedTable_;

#define makelimitedSurfaceInterpolationScheme(SS)                              \
                                                                               \
makelimitedSurfaceInterpolationTypeScheme(SS, scalar)                          \
makelimitedSurfaceInterpolationTypeScheme(SS, vector)                          \
makelimitedSurfaceInterpolationTypeScheme(SS, sphericalTensor)                 \
makelimitedSurfaceInterpolationTypeScheme(SS, symmTensor)                      \
makelimitedSurfaceInterpolationTypeScheme(SS, tensor)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "limitedSurfaceInterpolationScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
