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
    Foam::waveSpectrum

Description
    Base class for wave spectra

SourceFiles
    waveSpectrum.C

\*---------------------------------------------------------------------------*/

#ifndef waveSpectrum_H
#define waveSpectrum_H

#include "dictionary.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class waveSpectrum Declaration
\*---------------------------------------------------------------------------*/

class waveSpectrum
{
    // Private Data

        //- Gravitational acceleration [m/s^2]
        const scalar g_;

        //- Number of intervals to use for integration
        const label n_;


    // Private Member Functions

        //- Refine or coarsen the given field to the given size, using linear
        //  interpolation to construct intermediate values. Guaranteed
        //  symmetric, so refined(x.size(), refined(n, x)) will return x (if n
        //  is larger than x.size()).
        static tmp<scalarField> refined(const label n, const scalarField& x);

        //- Tmp variant of the above
        static tmp<scalarField> refined
        (
            const label n,
            const tmp<scalarField>& x
        );


protected:

    // Protected Member Functions

        //- Return the frequency below which a given fraction of the spectrum's
        //  total energy falls, given also the maximum frequency with
        //  appreciable energy
        scalar fFraction(const scalar fraction, const scalar f1) const;


public:

    //- Runtime type information
    TypeName("waveSpectrum");


    // Declare runtime construction
    declareRunTimeSelectionTable
    (
        autoPtr,
        waveSpectrum,
        dictionary,
        (const dictionary& dict, const scalar g),
        (dict, g)
    );


    // Constructors

        //- Construct a copy
        waveSpectrum(const waveSpectrum& spectrum);

        //- Construct from a dictionary and gravity
        waveSpectrum(const dictionary& dict, const scalar g);

        //- Construct a clone
        virtual autoPtr<waveSpectrum> clone() const = 0;


    // Selectors

        //- Select given a dictionary and gravity
        static autoPtr<waveSpectrum> New
        (
            const dictionary& dict,
            const scalar g
        );


    //- Destructor
    virtual ~waveSpectrum();


    // Member Functions

        //- Access the gravitation acceleration [m/s^2]
        inline scalar g() const
        {
            return g_;
        }

        //- Evaluate the wave spectral density at the given frequencies [m^2/Hz]
        virtual tmp<scalarField> S(const scalarField& f) const = 0;

        //- Evaluate the integral of the wave spectral density at the given
        //  frequencies [m^2]
        virtual tmp<scalarField> integralS(const scalarField& f) const;

        //- Evaluate the integral of the wave spectral density multiplied by
        //  the frequency at the given frequencies [Hz m^2]
        virtual tmp<scalarField> integralFS(const scalarField& f) const;

        //- Return the frequency below which a given fraction of the spectrum's
        //  total energy falls
        virtual scalar fFraction(const scalar fraction) const = 0;

        //- Write
        virtual void write(Ostream& os) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const waveSpectrum&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
