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
    Foam::CentredFitScheme

Description
    Centred fit surface interpolation scheme which applies an explicit
    correction to linear.

\*---------------------------------------------------------------------------*/

#ifndef CentredFitScheme_H
#define CentredFitScheme_H

#include "CentredFitData.H"
#include "linear.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class CentredFitScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class Polynomial, class Stencil>
class CentredFitScheme
:
    public linear<Type>
{
    // Private Data

        //- Factor the fit is allowed to deviate from linear.
        //  This limits the amount of high-order correction and increases
        //  stability on bad meshes
        const scalar linearLimitFactor_;

        //- Weights for central stencil
        const scalar centralWeight_;


public:

    //- Runtime type information
    TypeName("CentredFitScheme");


    // Constructors

        //- Construct from mesh and Istream
        CentredFitScheme(const fvMesh& mesh, Istream& is)
        :
            linear<Type>(mesh),
            linearLimitFactor_(readScalar(is)),
            centralWeight_(1000)
        {}


        //- Construct from mesh, faceFlux and Istream
        CentredFitScheme
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            linear<Type>(mesh),
            linearLimitFactor_(readScalar(is)),
            centralWeight_(1000)
        {}

        //- Disallow default bitwise copy construction
        CentredFitScheme(const CentredFitScheme&) = delete;


    // Member Functions

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return true;
        }

        //- Return the explicit correction to the face-interpolate
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        correction
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const
        {
            const fvMesh& mesh = this->mesh();

            const extendedCentredCellToFaceStencil& stencil = Stencil::New
            (
                mesh
            );

            const CentredFitData<Polynomial>& cfd =
            CentredFitData<Polynomial>::New
            (
                mesh,
                stencil,
                linearLimitFactor_,
                centralWeight_
            );

            const List<scalarList>& f = cfd.coeffs();

            return stencil.weightedSum(vf, f);
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const CentredFitScheme&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables

#define makeCentredFitSurfaceInterpolationTypeScheme\
(                                                                              \
    SS,                                                                        \
    POLYNOMIAL,                                                                \
    STENCIL,                                                                   \
    TYPE                                                                       \
)                                                                              \
                                                                               \
typedef CentredFitScheme<TYPE, POLYNOMIAL, STENCIL>                            \
    CentredFitScheme##TYPE##POLYNOMIAL##STENCIL##_;                            \
defineTemplateTypeNameAndDebugWithName                                         \
    (CentredFitScheme##TYPE##POLYNOMIAL##STENCIL##_, #SS, 0);                  \
                                                                               \
surfaceInterpolationScheme<TYPE>::addMeshConstructorToTable                    \
<CentredFitScheme<TYPE, POLYNOMIAL, STENCIL>>                                  \
    add##SS##STENCIL##TYPE##MeshConstructorToTable_;                           \
                                                                               \
surfaceInterpolationScheme<TYPE>::addMeshFluxConstructorToTable                \
<CentredFitScheme<TYPE, POLYNOMIAL, STENCIL>>                                  \
    add##SS##STENCIL##TYPE##MeshFluxConstructorToTable_;

#define makeCentredFitSurfaceInterpolationScheme(SS, POLYNOMIAL, STENCIL)      \
                                                                               \
makeCentredFitSurfaceInterpolationTypeScheme(SS,POLYNOMIAL,STENCIL,scalar)     \
makeCentredFitSurfaceInterpolationTypeScheme(SS,POLYNOMIAL,STENCIL,vector)     \
makeCentredFitSurfaceInterpolationTypeScheme                                   \
(                                                                              \
    SS,                                                                        \
    POLYNOMIAL,                                                                \
    STENCIL,                                                                   \
    sphericalTensor                                                            \
)                                                                              \
makeCentredFitSurfaceInterpolationTypeScheme(SS,POLYNOMIAL,STENCIL,symmTensor)\
makeCentredFitSurfaceInterpolationTypeScheme(SS,POLYNOMIAL,STENCIL,tensor)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
