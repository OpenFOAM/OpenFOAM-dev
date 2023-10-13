/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    Foam::waveModels::irregular

Description
    Irregular wave model. Phase fraction and velocity field are built up from
    multiple first-order waves, sampled from a selectable wave spectrum.

Usage
    \table
        Property  | Description                     | Required? | Default
        depth     | The water depth [m]             | no        | great
        spectrum  | The wave spectrum               | yes       |
        n         | The number of times to sample \\
                    the spectrum                    | yes       |
        span      | The fractional range across \\
                    which to sample the spectrum    | no        | (0.01 0.99)
        setFormat | The format with which to plot \\
                    the spectrum                    | no        | none
    \endtable

    Example specification in constant/waveProperties:
    \verbatim
    waves
    (
        irregular
        {
            spectrum    PiersonMoskowitz;

            PiersonMoskowitzCoeffs
            {
                U19_5       15;
            }

            n           16;
            span        (0.01 0.99);
        }
    );
    \endverbatim

SourceFiles
    irregular.C

\*---------------------------------------------------------------------------*/

#ifndef irregular_H
#define irregular_H

#include "waveModel.H"
#include "waveSpectrum.H"
#include "AiryCoeffs.H"
#include "setWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveModels
{

/*---------------------------------------------------------------------------*\
                            Class irregular Declaration
\*---------------------------------------------------------------------------*/

class irregular
:
    public waveModel
{
    // Private Data

        //- Depth [m]
        const scalar depth_;

        //- Spectrum
        autoPtr<waveSpectrum> spectrum_;

        //- The number of times to sample the spectrum
        const label n_;

        //- The fractional range across which to sample the spectrum
        const Pair<scalar> span_;

        //- Formatter with which to plot the spectrum
        autoPtr<setWriter> formatter_;

        //- Amplitudes [m]
        scalarList amplitudes_;

        //- Wavelengths [m]
        scalarList lengths_;

        //- Phase offsets [rad]
        scalarList phases_;


    // Private Member Function

        //- Return the i-th wave coefficients
        AiryCoeffs coeffs(const label i) const;


public:

    //- Runtime type information
    TypeName("irregular");


    // Constructors

        //- Construct a copy
        irregular(const irregular& wave);

        //- Construct from a dictionary and gravity
        irregular
        (
            const dictionary& dict,
            const scalar g,
            const word& modelName = irregular::typeName
        );

        //- Construct a clone
        virtual autoPtr<waveModel> clone() const
        {
            return autoPtr<waveModel>(new irregular(*this));
        }


    //- Destructor
    virtual ~irregular();


    // Member Functions

        // Access

            //- Get the depth
            scalar depth() const
            {
                return depth_;
            }


        //- The wave celerity [m/s]
        virtual scalar celerity() const;

        //- Get the wave elevation at a given time and local coordinates. Local
        //  x is aligned with the direction of propagation.
        virtual tmp<scalarField> elevation
        (
            const scalar t,
            const scalarField& x
        ) const;

        //- Get the wave velocity at a given time and local coordinates. Local
        //  x is aligned with the direction of propagation, and z with negative
        //  gravity.
        virtual tmp<vector2DField> velocity
        (
            const scalar t,
            const vector2DField& xz
        ) const;

        //- Write
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
