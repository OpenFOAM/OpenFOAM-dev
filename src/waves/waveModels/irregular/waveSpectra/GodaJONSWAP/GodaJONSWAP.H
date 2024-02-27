/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    Foam::GodaJONSWAP

Description
    GodaJONSWAP wave spectrum. This is an alternative, approximate
    parameterisation of the JONSWAP spectrum, in which the significant wave
    height and period are specified instead of the wind speed and fetch.

    References:
    \verbatim
        Goda, Y. (1988).
        Statistical variability of sea state parameters as a function of wave
        spectrum.
        Coastal Engineering in Japan, 31(1), 39-52.
    \endverbatim

    \verbatim
        Goda, Y. (2010).
        Random seas and design of maritime structures.
        World Scientific Publishing Company.
    \endverbatim

    See page 29 of the second reference for a convenient formulation.

Usage
    \table
        Property | Description                             | Required? | Default
        Hs       | The significant wave height [m]         | yes       |
        Tp       | The significant wave period [s]         | yes       |
        gamma    | Peaked-ness parameter                   | no        | 3.3
    \endtable

    Example specification:
    \verbatim
    spectrum    GodaJONSWAP;

    GodaJONSWAPCoeffs
    {
        Hs          2;
        Ts          6;
    }
    \endverbatim

See also
    Foam::waveSpectra::JONSWAP

SourceFiles
    GodaJONSWAP.C

\*---------------------------------------------------------------------------*/

#ifndef GodaJONSWAP_H
#define GodaJONSWAP_H

#include "waveSpectrum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveSpectra
{

/*---------------------------------------------------------------------------*\
                          Class GodaJONSWAP Declaration
\*---------------------------------------------------------------------------*/

class GodaJONSWAP
:
    public waveSpectrum
{
    // Private Data

        //- Significant Wave Height [m]
        const scalar Hs_;

        //- Significant wave period [s]
        const scalar Ts_;

        //- Peaked-ness parameter
        const scalar gamma_;


public:

    //- Runtime type information
    TypeName("GodaJONSWAP");


    // Constructors

        //- Construct a copy
        GodaJONSWAP(const GodaJONSWAP& spectrum);

        //- Construct from a dictionary and gravity
        GodaJONSWAP(const dictionary& dict, const scalar g);

        //- Construct a clone
        virtual autoPtr<waveSpectrum> clone() const
        {
            return autoPtr<waveSpectrum>(new GodaJONSWAP(*this));
        }


    //- Destructor
    virtual ~GodaJONSWAP();


    // Member Functions

        //- Evaluate the wave spectral density at the given frequencies [m^2/Hz]
        virtual tmp<scalarField> S(const scalarField& f) const;

        //- Return the frequency below which a given fraction of the spectrum's
        //  total energy falls []
        virtual scalar fFraction(const scalar fraction) const;

        //- Write
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveSpectra
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
