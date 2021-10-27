/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2021 OpenFOAM Foundation
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
    Foam::interfaceCompression

Description
    Interface compression corrected scheme, based on counter-gradient
    transport, to maintain sharp interfaces during VoF simulations.

    The interface compression is applied to the face interpolated field from a
    suitable 2nd-order shape-preserving NVD or TVD scheme, e.g.  vanLeer or
    vanAlbada.  A coefficient is supplied to control the degree of compression,
    with a value of 1 suitable for most VoF cases to ensure interface integrity.
    A value larger than 1 can be used but the additional compression can bias
    the interface to follow the mesh more closely while a value smaller than 1
    can lead to interface smearing.

    Example:
    \verbatim
    divSchemes
    {
        .
        .
        div(phi,alpha)     Gauss interfaceCompression vanLeer 1;
        .
        .
    }
    \endverbatim

See also
    Foam::vanLeer
    Foam::vanAlbada
    Foam::PLIC
    Foam::PLICU
    Foam::MPLIC
    Foam::MPLICU

SourceFiles
    interfaceCompression.C

\*---------------------------------------------------------------------------*/

#ifndef interfaceCompression_H
#define interfaceCompression_H

#include "surfaceInterpolationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class interfaceCompression Declaration
\*---------------------------------------------------------------------------*/

class interfaceCompressionNew
:
    public surfaceInterpolationScheme<scalar>
{
    // Private member data

        const surfaceScalarField& phi_;

        //- Base scheme to which the compression is applied
        tmp<surfaceInterpolationScheme<scalar>> tScheme_;

        //- Compression factor
        const scalar cAlpha_;


public:

    //- Runtime type information
    TypeName("interfaceCompression");


    // Constructors

        //- Construct from faceFlux and Istream
        interfaceCompressionNew
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            surfaceInterpolationScheme<scalar>(mesh),
            phi_(faceFlux),
            tScheme_
            (
                surfaceInterpolationScheme<scalar>::New(mesh, faceFlux, is)
            ),
            cAlpha_(readScalar(is))
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
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const interfaceCompressionNew&) = delete;
};


extern const wordHashSet compressionSchemes;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
