/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::noiseFFT

Description
    FFT of the pressure field

SourceFiles
    noiseFFT.C

\*---------------------------------------------------------------------------*/

#ifndef noiseFFT_H
#define noiseFFT_H

#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class noiseFFT Declaration
\*---------------------------------------------------------------------------*/

class noiseFFT
:
    public scalarField
{
    // Private Data

        //- Time spacing of the raw data
        scalar deltat_;


public:

    //- Reference pressure
    static scalar p0;


    // Constructors

        //- Construct from pressure field
        noiseFFT
        (
            const scalar deltat,
            const scalarField& pressure
        );

        //- Construct from Istream
        noiseFFT(Istream&);

        //- Construct from pressure field file name
        noiseFFT(const fileName& pFileName, const label skip = 0);


    // Member Functions

        //- Return p(t)
        Pair<scalarField> pt() const;

        //- Return the nth window
        tmp<scalarField> window(const label N, const label n) const;

        //- Return the Hanning window function
        tmp<scalarField> Hanning(const label N) const;

        //- Return the fft of the given pressure data
        tmp<scalarField> Pf(const tmp<scalarField>& pn) const;

        //- Return the multi-window mean fft of the complete pressure data
        Pair<scalarField> meanPf(const label N, const label nw) const;

        //- Return the multi-window RMS mean fft of the complete pressure data
        Pair<scalarField> RMSmeanPf(const label N, const label nw) const;

        //- Return the narrow-band PFL (pressure-fluctuation level) spectrum
        Pair<scalarField> Lf(const Pair<scalarField>& gPf) const;

        //- Return the one-third-octave-band PFL spectrum
        //  starting at octave with mean frequency f1
        Pair<scalarField> Ldelta
        (
            const Pair<scalarField>& gLf,
            const scalar f1,
            const scalar fU
        ) const;

        //- Return the one-third-octave-band pressure spectrum
        //  starting at octave with mean frequency f1
        Pair<scalarField> Pdelta
        (
            const Pair<scalarField>& gLf,
            const scalar f1,
            const scalar fU
        ) const;

        //- Return the total PFL as the sum of Lf over all frequencies
        scalar Lsum(const Pair<scalarField>& gLf) const;

        //- Convert the db into Pa
        scalar dbToPa(const scalar db) const;

        //- Convert the db-field into Pa
        tmp<scalarField> dbToPa(const tmp<scalarField>& db) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
