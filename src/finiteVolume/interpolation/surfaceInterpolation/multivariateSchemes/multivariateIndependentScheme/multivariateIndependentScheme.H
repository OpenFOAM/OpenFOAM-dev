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
    Foam::multivariateIndependentScheme

Description
    Generic multi-variate discretisation scheme class for which any of the
    NVD, CNVD or NVDV schemes may be selected for each variable and applied
    independently.

    This is equivalent to using separate "div" terms and schemes for each
    variable/equation.

SourceFiles
    multivariateIndependentScheme.C

\*---------------------------------------------------------------------------*/

#ifndef multivariateIndependentScheme_H
#define multivariateIndependentScheme_H

#include "multivariateSurfaceInterpolationScheme.H"
#include "limitedSurfaceInterpolationScheme.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class multivariateIndependentScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class multivariateIndependentScheme
:
    public multivariateSurfaceInterpolationScheme<Type>
{
    // Private data

        dictionary schemes_;
        const surfaceScalarField& faceFlux_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        multivariateIndependentScheme(const multivariateIndependentScheme&);

        //- Disallow default bitwise assignment
        void operator=(const multivariateIndependentScheme&);


public:

    //- Runtime type information
    TypeName("multivariateIndependent");


    // Constructors

        //- Construct for field, faceFlux and Istream
        multivariateIndependentScheme
        (
            const fvMesh& mesh,
            const typename multivariateSurfaceInterpolationScheme<Type>::
                fieldTable& fields,
            const surfaceScalarField& faceFlux,
            Istream& schemeData
        );


    // Member Operators

        tmp<surfaceInterpolationScheme<Type>> operator()
        (
            const GeometricField<Type, fvPatchField, volMesh>& field
        ) const
        {
            return surfaceInterpolationScheme<Type>::New
            (
                faceFlux_.mesh(),
                faceFlux_,
                schemes_.lookup(field.name())
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "multivariateIndependentScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
