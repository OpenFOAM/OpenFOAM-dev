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
    Foam::PiersonMoskowitz

Description
    Pierson-Moskowitz wave spectrum. This spectrum has the following form:

        \f[
            S(\omega) = \frac{\alpha g^2}{\omega^5} \exp \left(- \beta \\
            \left( \frac{\omega_0}{\omega} \right)^4 \right)
        \f]

    \vartable
        \omega    | angular frequency [rad/s], equal to \f$2 \pi f\f$
        f         | frequency [Hz]
        S(\omega) | spectral density [m^2/Hz]
        \alpha    | coefficient, equal to 8.1e3
        \beta     | coefficient, equal to 0.74
        \omega_0  | reference angular frequency, equal to \f$g/U_{19.5}\f$
        U_{19.5}  | reference velocity magnitude, 19.5 metres above the sea \\
                    surface
    \endvartable

    References:
    \verbatim
        Pierson Jr, W. J., & Moskowitz, L. (1964).
        A proposed spectral form for fully developed wind seas based on the \\
        similarity theory of SA Kitaigorodskii.
        Journal of geophysical research, 69(24), 5181-5190.
    \endverbatim

    \verbatim
        Stewart, R. H. (2008).
        Introduction to physical oceanography.
        Robert H. Stewart.
    \endverbatim

Usage
    \table
        Property | Description                             | Required? | Default
        U19_5    | The air speed 19.5 metres above the \\
                   surface [m/s]                           | yes       |
    \endtable

    Example specification:
    \verbatim
    spectrum    PiersonMoskowitz;

    PiersonMoskowitzCoeffs
    {
        U19_5       15;
    }
    \endverbatim

SourceFiles
    PiersonMoskowitz.C

\*---------------------------------------------------------------------------*/

#ifndef PiersonMoskowitz_H
#define PiersonMoskowitz_H

#include "waveSpectrum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveSpectra
{

/*---------------------------------------------------------------------------*\
                          Class PiersonMoskowitz Declaration
\*---------------------------------------------------------------------------*/

class PiersonMoskowitz
:
    public waveSpectrum
{
    // Private Data

        //- Reference velocity magnitude
        const scalar U19_5_;

        //- Alpha coefficient
        const scalar alpha_;

        //- Beta coefficient
        const scalar beta_;


public:

    //- Runtime type information
    TypeName("PiersonMoskowitz");


    // Constructors

        //- Construct a copy
        PiersonMoskowitz(const PiersonMoskowitz& spectrum);

        //- Construct from a dictionary and gravity
        PiersonMoskowitz(const dictionary& dict, const scalar g);

        //- Construct a clone
        virtual autoPtr<waveSpectrum> clone() const
        {
            return autoPtr<waveSpectrum>(new PiersonMoskowitz(*this));
        }


    //- Destructor
    virtual ~PiersonMoskowitz();


    // Member Functions

        //- Evaluate the wave spectral density at the given frequencies [m^2/Hz]
        virtual tmp<scalarField> S(const scalarField& f) const;

        //- Evaluate the integral of the wave spectral density at the given
        //  frequencies [m^2]
        virtual tmp<scalarField> integralS(const scalarField& f) const;

        //- Return the frequency below which a given fraction of the spectrum's
        //  total power falls []
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
