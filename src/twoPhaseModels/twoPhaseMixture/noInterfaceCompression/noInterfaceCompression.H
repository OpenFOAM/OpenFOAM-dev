/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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
    Foam::noInterfaceCompression

Description
    Wrapper scheme to allow VoF solvers to run efficiently without interface
    compression, e.g. for cavitation simulations.

    Example:
    \verbatim
    divSchemes
    {
        .
        .
        div(phi,alpha)     Gauss noInterfaceCompression vanLeer;
        .
        .
    }
    \endverbatim

See also
    Foam::interfaceCompression
    Foam::vanLeer
    Foam::vanAlbada

SourceFiles
    noInterfaceCompression.C

\*---------------------------------------------------------------------------*/

#ifndef noInterfaceCompression_H
#define noInterfaceCompression_H

#include "surfaceInterpolationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class noInterfaceCompression Declaration
\*---------------------------------------------------------------------------*/

class noInterfaceCompressionNew
:
    public surfaceInterpolationScheme<scalar>
{
    // Private member data

        //- Base scheme to which the compression is applied
        tmp<surfaceInterpolationScheme<scalar>> tScheme_;


public:

    //- Runtime type information
    TypeName("noInterfaceCompression");


    // Constructors

        //- Construct from faceFlux and Istream
        noInterfaceCompressionNew
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            surfaceInterpolationScheme<scalar>(mesh),
            tScheme_
            (
                surfaceInterpolationScheme<scalar>::New(mesh, faceFlux, is)
            )
        {}


    // Member Functions

        //- Return the interpolation weighting factors
        virtual tmp<surfaceScalarField> weights
        (
            const volScalarField&
        ) const
        {
            NotImplemented;

            return tmp<surfaceScalarField>(nullptr);
        }

        //- Return the face-interpolate of the given cell field
        virtual tmp<surfaceScalarField> interpolate
        (
            const volScalarField& vf
        ) const
        {
            return tScheme_().interpolate(vf);
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const noInterfaceCompressionNew&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
