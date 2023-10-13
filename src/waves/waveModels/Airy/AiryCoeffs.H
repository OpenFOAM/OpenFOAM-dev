/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2023 OpenFOAM Foundation
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
    Foam::waveModels::AiryCoeffs

Description
    Calculation engine for the Airy wave model and other models that are a
    correction on top of the Airy model or a superposition of Airy models

SourceFiles
    AiryCoeffs.C

\*---------------------------------------------------------------------------*/

#ifndef AiryCoeffs_H
#define AiryCoeffs_H

#include "waveModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveModels
{

/*---------------------------------------------------------------------------*\
                           Class AiryCoeffs Declaration
\*---------------------------------------------------------------------------*/

class AiryCoeffs
{
private:

    // Private Static Member Functions

        //- Calculate the length from the depth, amplitude, period and g.
        //  Requires a function which returns celerity given the depth,
        //  amplitude, length and g. Iterates to solve this function for the
        //  specified period.
        static scalar calcLength
        (
            const scalar depth,
            const scalar amplitude,
            const scalar period,
            const scalar g,
            scalar (*celerityPtr)(const AiryCoeffs&)
        );


public:

    // Public Data

        //- Depth [m]
        const scalar depth;

        //- Amplitude [m]
        const scalar amplitude;

        //- Wavelength [m]
        const scalar length;

        //- Gravitational acceleration [m/s^2]
        const scalar g;


    // Constructors

        //- Construct from components
        AiryCoeffs
        (
            const scalar depth,
            const scalar amplitude,
            const scalar length,
            const scalar g
        );

        //- Construct from components but with period instead of length
        AiryCoeffs
        (
            const scalar depth,
            const scalar amplitude,
            const scalar period,
            const scalar g,
            scalar (*celerityPtr)(const AiryCoeffs&)
        );


    // Member Functions

        //- The angular wavenumber [rad/m]
        scalar k() const;

        //- Return whether shallow and intermediate effects are to be omitted
        bool deep() const;

        //- The wave celerity [m/s]
        static scalar celerity(const AiryCoeffs& coeffs);

        //- The wave celerity [m/s]
        scalar celerity() const;

        //- Angle of the oscillation [rad]
        tmp<scalarField> angle
        (
            const scalar phase,
            const scalar t,
            const scalarField& x
        ) const;

        //- Get the wave elevation at a given time and local coordinates. Local
        //  x is aligned with the direction of propagation.
        tmp<scalarField> elevation
        (
            const scalar phase,
            const scalar t,
            const scalarField& x
        ) const;

        //- Return the non-dimensionalised i-th harmonic of the velocity
        tmp<vector2DField> vi
        (
            const label i,
            const scalar phase,
            const scalar t,
            const vector2DField& xz
        ) const;

        //- Get the wave velocity at a given time and local coordinates. Local
        //  x is aligned with the direction of propagation, and z with negative
        //  gravity.
        tmp<vector2DField> velocity
        (
            const scalar phase,
            const scalar t,
            const vector2DField& xz
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
