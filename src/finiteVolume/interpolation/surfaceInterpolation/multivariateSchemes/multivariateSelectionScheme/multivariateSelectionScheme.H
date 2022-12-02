/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    Foam::multivariateSelectionScheme

Description
    Generic multi-variate discretisation scheme class for which any of the
    NVD, CNVD or NVDV schemes may be selected for each variable.

SourceFiles
    multivariateSelectionScheme.C

\*---------------------------------------------------------------------------*/

#ifndef multivariateSelectionScheme_H
#define multivariateSelectionScheme_H

#include "multivariateSurfaceInterpolationScheme.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class multivariateSelectionScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class multivariateSelectionScheme
:
    public multivariateSurfaceInterpolationScheme<Type>
{
    // Private Data

        dictionary schemes_;
        const surfaceScalarField& faceFlux_;
        surfaceScalarField weights_;


public:

    //- Runtime type information
    TypeName("multivariateSelection");


    // Constructors

        //- Construct for field, faceFlux and Istream
        multivariateSelectionScheme
        (
            const fvMesh& mesh,
            const typename multivariateSurfaceInterpolationScheme<Type>::
                fieldTable& fields,
            const surfaceScalarField& faceFlux,
            Istream& schemeData
        );

        //- Disallow default bitwise copy construction
        multivariateSelectionScheme
        (
            const multivariateSelectionScheme&
        ) = delete;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const multivariateSelectionScheme&) = delete;

        //- surfaceInterpolationScheme sub-class returned by operator(field)
        //  which is used as the interpolation scheme for the field
        class fieldScheme
        :
            public multivariateSurfaceInterpolationScheme<Type>::fieldScheme
        {
            // Private Data

                const surfaceScalarField& weights_;

        public:

            // Constructors

                fieldScheme
                (
                    const VolField<Type>& field,
                    const surfaceScalarField& weights
                )
                :
                    multivariateSurfaceInterpolationScheme<Type>::
                        fieldScheme(field),
                    weights_(weights)
                {}


            // Member Functions

                //- Return the interpolation weighting factors
                tmp<surfaceScalarField> weights
                (
                    const VolField<Type>&
                ) const
                {
                    return weights_;
                }
        };

        tmp<surfaceInterpolationScheme<Type>> operator()
        (
            const VolField<Type>& field
        ) const
        {
            return tmp<surfaceInterpolationScheme<Type>>
            (
                new fieldScheme(field, weights_)
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "multivariateSelectionScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
