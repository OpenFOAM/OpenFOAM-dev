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
    Foam::JONSWAP

Description
    JONSWAP wave spectrum. This is similar to the Pierson-Moskowitz spectrum,
    but with an additional empirical correction to account for the fetch
    (distance to the lee shore).

    References:
    \verbatim
        Hasselmann, K., Barnett, T. P., Bouws, E., Carlson, H., Cartwright, \\
        D. E., Enke, K., ... & Walden, H. (1973).
        Measurements of wind-wave growth and swell decay during the Joint \\
        North Sea Wave Project (JONSWAP).
        Ergaenzungsheft zur Deutschen Hydrographischen Zeitschrift, Reihe A.
    \endverbatim

    \verbatim
        Stewart, R. H. (2008).
        Introduction to physical oceanography.
        Robert H. Stewart.
    \endverbatim

Usage
    \table
        Property | Description                             | Required? | Default
        U10      | The air speed 10 metres above the \\
                   surface [m/s]                           | yes       |
        F        | The fetch (distance from the lee \\
                   shore) [m]                              | yes       |
        gamma    | Peaked-ness parameter                   | no        | 3.3
    \endtable

    Example specification:
    \verbatim
    spectrum    JONSWAP;

    JONSWAPCoeffs
    {
        U10         10;
        F           200e3;
    }
    \endverbatim

See also
    Foam::waveSpectra::PiersonMoskowitz

SourceFiles
    JONSWAP.C

\*---------------------------------------------------------------------------*/

#ifndef JONSWAP_H
#define JONSWAP_H

#include "waveSpectrum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveSpectra
{

/*---------------------------------------------------------------------------*\
                          Class JONSWAP Declaration
\*---------------------------------------------------------------------------*/

class JONSWAP
:
    public waveSpectrum
{
    // Private Data

        //- Reference velocity magnitude [m/s]
        const scalar U10_;

        //- Fetch (distance from the lee shore) [m]
        const scalar F_;

        //- Peaked-ness parameter
        const scalar gamma_;


public:

    //- Runtime type information
    TypeName("JONSWAP");


    // Constructors

        //- Construct a copy
        JONSWAP(const JONSWAP& spectrum);

        //- Construct from a dictionary and gravity
        JONSWAP(const dictionary& dict, const scalar g);

        //- Construct a clone
        virtual autoPtr<waveSpectrum> clone() const
        {
            return autoPtr<waveSpectrum>(new JONSWAP(*this));
        }


    //- Destructor
    virtual ~JONSWAP();


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
