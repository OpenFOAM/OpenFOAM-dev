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
    Foam::harmonic

Description
    Harmonic-mean interpolation scheme class.

    This scheme interpolates 1/field using a scheme specified at run-time
    and return the reciprocal of the interpolate.

SourceFiles
    harmonic.C

\*---------------------------------------------------------------------------*/

#ifndef harmonic_H
#define harmonic_H

#include "surfaceInterpolationScheme.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "reverseLinear.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class harmonic Declaration
\*---------------------------------------------------------------------------*/

class harmonic
:
    public surfaceInterpolationScheme<scalar>
{

public:

    //- Runtime type information
    TypeName("harmonic");


    // Constructors

        //- Construct from mesh
        harmonic(const fvMesh& mesh)
        :
            surfaceInterpolationScheme<scalar>(mesh)
        {}

        //- Construct from Istream.
        //  The name of the flux field is read from the Istream and looked-up
        //  from the mesh objectRegistry
        harmonic
        (
            const fvMesh& mesh,
            Istream& is
        )
        :
            surfaceInterpolationScheme<scalar>(mesh)
        {}

        //- Construct from faceFlux and Istream
        harmonic
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            surfaceInterpolationScheme<scalar>(mesh)
        {}


    // Member Functions

        //- Return the interpolation weighting factors
        virtual tmp<surfaceScalarField> weights
        (
            const GeometricField<scalar, fvPatchField, volMesh>&
        ) const
        {
            NotImplemented;

            return tmp<surfaceScalarField>(nullptr);
        }

        //- Return the face-interpolate of the given cell field
        virtual tmp<GeometricField<scalar, fvsPatchField, surfaceMesh>>
        interpolate
        (
            const GeometricField<scalar, fvPatchField, volMesh>& vf
        ) const
        {
            return 1.0/(reverseLinear<scalar>(vf.mesh()).interpolate(1.0/vf));
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const harmonic&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
