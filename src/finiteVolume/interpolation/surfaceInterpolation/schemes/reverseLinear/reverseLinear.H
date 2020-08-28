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
    Foam::reverseLinear

Description
    Inversed weight central-differencing interpolation scheme class.

    Useful for inverse weighted and harmonic interpolations.

SourceFiles
    reverseLinear.C

\*---------------------------------------------------------------------------*/

#ifndef reverseLinear_H
#define reverseLinear_H

#include "surfaceInterpolationScheme.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class reverseLinear Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class reverseLinear
:
    public surfaceInterpolationScheme<Type>
{
public:

    //- Runtime type information
    TypeName("reverseLinear");


    // Constructors

        //- Construct from mesh
        reverseLinear(const fvMesh& mesh)
        :
            surfaceInterpolationScheme<Type>(mesh)
        {}

        //- Construct from Istream
        reverseLinear(const fvMesh& mesh, Istream&)
        :
            surfaceInterpolationScheme<Type>(mesh)
        {}

        //- Construct from faceFlux and Istream
        reverseLinear
        (
            const fvMesh& mesh,
            const surfaceScalarField&,
            Istream&
        )
        :
            surfaceInterpolationScheme<Type>(mesh)
        {}


    // Member Functions

        //- Return the interpolation weighting factors
        tmp<surfaceScalarField> weights
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const
        {
            const fvMesh& mesh = this->mesh();

            tmp<surfaceScalarField> tcdWeights
            (
                mesh.surfaceInterpolation::weights()
            );
            const surfaceScalarField& cdWeights = tcdWeights();

            tmp<surfaceScalarField> treverseLinearWeights
            (
                surfaceScalarField::New
                (
                    "reverseLinearWeights",
                    mesh,
                    dimless
                )
            );
            surfaceScalarField& reverseLinearWeights =
                treverseLinearWeights.ref();

            reverseLinearWeights.primitiveFieldRef() =
                1.0 - cdWeights.primitiveField();

            surfaceScalarField::Boundary& rlwbf =
                reverseLinearWeights.boundaryFieldRef();


            forAll(mesh.boundary(), patchi)
            {
                if (rlwbf[patchi].coupled())
                {
                    rlwbf[patchi] = 1.0 - cdWeights.boundaryField()[patchi];
                }
                else
                {
                    rlwbf[patchi] = cdWeights.boundaryField()[patchi];
                }
            }

            return treverseLinearWeights;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const reverseLinear&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
