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
    Foam::fv::faceMDLimitedGrad

Description
    faceMDLimitedGrad gradient scheme applied to a runTime selected
    base gradient scheme.

    The scalar limiter based on limiting the extrapolated face values
    between the face-neighbour cell values and is applied to the gradient
    in each face direction separately.

SourceFiles
    faceMDLimitedGrad.C

\*---------------------------------------------------------------------------*/

#ifndef faceMDLimitedGrad_H
#define faceMDLimitedGrad_H

#include "gradScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                      Class faceMDLimitedGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class faceMDLimitedGrad
:
    public fv::gradScheme<Type>
{
    // Private Data

        tmp<fv::gradScheme<Type>> basicGradScheme_;

        //- Limiter coefficient
        const scalar k_;


    // Private Member Functions

        inline void limitFace
        (
            scalar& limiter,
            const scalar maxDelta,
            const scalar minDelta,
            const scalar extrapolate
        ) const;


        //- Disallow default bitwise copy construction
        faceMDLimitedGrad(const faceMDLimitedGrad&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const faceMDLimitedGrad&) = delete;


public:

    //- RunTime type information
    TypeName("faceMDLimited");


    // Constructors

        //- Construct from mesh and schemeData
        faceMDLimitedGrad(const fvMesh& mesh, Istream& schemeData)
        :
            gradScheme<Type>(mesh),
            basicGradScheme_(fv::gradScheme<Type>::New(mesh, schemeData)),
            k_(readScalar(schemeData))
        {
            if (k_ < 0 || k_ > 1)
            {
                FatalIOErrorInFunction
                (
                    schemeData
                )   << "coefficient = " << k_
                    << " should be >= 0 and <= 1"
                    << exit(FatalIOError);
            }
        }


    // Member Functions

        //- Return the gradient of the given field to the gradScheme::grad
        //  for optional caching
        virtual tmp
        <
            GeometricField
            <typename outerProduct<vector, Type>::type, fvPatchField, volMesh>
        > calcGrad
        (
            const GeometricField<Type, fvPatchField, volMesh>& vsf,
            const word& name
        ) const;
};


// * * * * * * * * Template Member Function Specialisations  * * * * * * * * //

template<>
tmp<volVectorField> faceMDLimitedGrad<scalar>::calcGrad
(
    const volScalarField& vsf,
    const word& name
) const;


template<>
tmp<volTensorField> faceMDLimitedGrad<vector>::calcGrad
(
    const volVectorField& vsf,
    const word& name
) const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
