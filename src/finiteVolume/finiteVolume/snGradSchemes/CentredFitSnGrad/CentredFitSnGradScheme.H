/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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
    Foam::fv::CentredFitSnGradScheme

Description
    Centred fit snGrad scheme which applies an explicit correction to snGrad

\*---------------------------------------------------------------------------*/

#ifndef CentredFitSnGradScheme_H
#define CentredFitSnGradScheme_H

#include "CentredFitSnGradData.H"
#include "snGradScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{
/*---------------------------------------------------------------------------*\
                           Class CentredFitSnGradSnGradScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class Polynomial, class Stencil>
class CentredFitSnGradScheme
:
    public snGradScheme<Type>
{
    // Private Data

        //- Factor the fit is allowed to deviate from linear.
        //  This limits the amount of high-order correction and increases
        //  stability on bad meshes
        const scalar linearLimitFactor_;

        //- Weights for central stencil
        const scalar centralWeight_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        CentredFitSnGradScheme(const CentredFitSnGradScheme&);

        //- Disallow default bitwise assignment
        void operator=(const CentredFitSnGradScheme&);


public:

    //- Runtime type information
    TypeName("CentredFitSnGradScheme");


    // Constructors

        //- Construct from mesh and Istream
        CentredFitSnGradScheme(const fvMesh& mesh, Istream& is)
        :
            snGradScheme<Type>(mesh),
            linearLimitFactor_(readScalar(is)),
            centralWeight_(1000)
        {}


    // Member Functions

        //- Return the interpolation weighting factors for the given field
        virtual tmp<surfaceScalarField> deltaCoeffs
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const
        {
            return this->mesh().nonOrthDeltaCoeffs();
        }

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

            const CentredFitSnGradData<Polynomial>& cfd =
                CentredFitSnGradData<Polynomial>::New
                (
                    mesh,
                    stencil,
                    linearLimitFactor_,
                    centralWeight_
                );

            tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> sft
            (
                stencil.weightedSum(vf, cfd.coeffs())
            );

            sft.ref().dimensions() /= dimLength;

            return sft;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables

#define makeCentredFitSnGradTypeScheme(SS, POLYNOMIAL, STENCIL, TYPE)          \
    typedef Foam::fv::CentredFitSnGradScheme                                   \
        <Foam::TYPE, Foam::POLYNOMIAL, Foam::STENCIL>                          \
        CentredFitSnGradScheme##TYPE##POLYNOMIAL##STENCIL##_;                  \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
        (CentredFitSnGradScheme##TYPE##POLYNOMIAL##STENCIL##_, #SS, 0);        \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        namespace fv                                                           \
        {                                                                      \
            snGradScheme<TYPE>::addMeshConstructorToTable                      \
                <CentredFitSnGradScheme<TYPE, POLYNOMIAL, STENCIL>>            \
                add##SS##STENCIL##TYPE##MeshConstructorToTable_;               \
        }                                                                      \
    }

#define makeCentredFitSnGradScheme(SS, POLYNOMIAL, STENCIL)                    \
                                                                               \
    makeCentredFitSnGradTypeScheme(SS,POLYNOMIAL,STENCIL,scalar)               \
    makeCentredFitSnGradTypeScheme(SS,POLYNOMIAL,STENCIL,vector)               \
    makeCentredFitSnGradTypeScheme(SS,POLYNOMIAL,STENCIL,sphericalTensor)      \
    makeCentredFitSnGradTypeScheme(SS,POLYNOMIAL,STENCIL,symmTensor)           \
    makeCentredFitSnGradTypeScheme(SS,POLYNOMIAL,STENCIL,tensor)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
