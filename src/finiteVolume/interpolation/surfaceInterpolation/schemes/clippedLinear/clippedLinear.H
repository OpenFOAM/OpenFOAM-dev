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
    Foam::clippedLinear

Description
    Central-differencing interpolation scheme using clipped-weights to
    improve stability on meshes with very rapid variations in cell size.

SourceFiles
    clippedLinear.C

\*---------------------------------------------------------------------------*/

#ifndef clippedLinear_H
#define clippedLinear_H

#include "surfaceInterpolationScheme.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class clippedLinear Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class clippedLinear
:
    public surfaceInterpolationScheme<Type>
{
    // Private data

        const scalar cellSizeRatio_;
        scalar wfLimit_;


    // Private Member Functions

        void calcWfLimit()
        {
            if (cellSizeRatio_ <= 0 || cellSizeRatio_ > 1)
            {
                FatalErrorInFunction
                    << "Given cellSizeRatio of " << cellSizeRatio_
                    << " is not between 0 and 1"
                    << exit(FatalError);
            }

            wfLimit_ = cellSizeRatio_/(1.0 + cellSizeRatio_);
        }


        //- Disallow default bitwise assignment
        void operator=(const clippedLinear&);


public:

    //- Runtime type information
    TypeName("clippedLinear");


    // Constructors

        //- Construct from mesh and cellSizeRatio
        clippedLinear(const fvMesh& mesh, const scalar cellSizeRatio)
        :
            surfaceInterpolationScheme<Type>(mesh),
            cellSizeRatio_(cellSizeRatio)
        {
            calcWfLimit();
        }

        //- Construct from Istream
        clippedLinear(const fvMesh& mesh, Istream& is)
        :
            surfaceInterpolationScheme<Type>(mesh),
            cellSizeRatio_(readScalar(is))
        {
            calcWfLimit();
        }

        //- Construct from faceFlux and Istream
        clippedLinear
        (
            const fvMesh& mesh,
            const surfaceScalarField&,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            cellSizeRatio_(readScalar(is))
        {
            calcWfLimit();
        }


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

            tmp<surfaceScalarField> tclippedLinearWeights
            (
                new surfaceScalarField
                (
                    IOobject
                    (
                        "clippedLinearWeights",
                        mesh.time().timeName(),
                        mesh
                    ),
                    mesh,
                    dimless
                )
            );
            surfaceScalarField& clippedLinearWeights =
                tclippedLinearWeights.ref();

            clippedLinearWeights.primitiveFieldRef() =
                max(min(cdWeights.primitiveField(), 1 - wfLimit_), wfLimit_);

            surfaceScalarField::Boundary& clwbf =
                clippedLinearWeights.boundaryFieldRef();

            forAll(mesh.boundary(), patchi)
            {
                if (clwbf[patchi].coupled())
                {
                    clwbf[patchi] =
                        max
                        (
                            min
                            (
                                cdWeights.boundaryField()[patchi],
                                1 - wfLimit_
                            ),
                            wfLimit_
                        );
                }
                else
                {
                    clwbf[patchi] = cdWeights.boundaryField()[patchi];
                }
            }

            return tclippedLinearWeights;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
