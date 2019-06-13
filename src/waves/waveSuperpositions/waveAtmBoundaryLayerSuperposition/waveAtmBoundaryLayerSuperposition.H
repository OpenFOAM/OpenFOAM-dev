/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2019 OpenFOAM Foundation
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
    Foam::waveAtmBoundaryLayerSuperposition

Description
    An extension of waveSuperposition which adds an atmospheric boundary layer
    model to the gas velocity. The user supplies a gas velocity and a height
    above the wave coordinate system origin at which that velocity is reached.
    Also needed are a maximum and minimum wave height which are used to set the
    surface roughness in the boundary layer model. It is not trivial to
    determine these from an arbitrary superposition of differently oriented
    wave models, so they are required as user inputs instead. For a pure
    sinusoidal wave, the maximum and minimum wave heights can be set to
    positive and negative amplitude, respectively.

Usage
    \table
        Property | Description                              | Req'd? | Default
        UGasRef  | The gas velocity at the reference height | yes    |
        hRef     | The reference height relative to the \\
                   origin of the wave coordinate system     | yes    |
        hWaveMin | The minimum wave elevation               | yes    |
        hWaveMax | The maximum wave elevation               | yes    |
    \endtable

    Example specification:
    \verbatim
    type        waveAtmBoundaryLayer;
    // parameters for waveSuperposition ...
    UGasRef     (10 0 0);
    hRef        10;
    hWaveMin    -2;
    hWaveMax    3;
    \endverbatim

SourceFiles
    waveAtmBoundaryLayerSuperposition.C

\*---------------------------------------------------------------------------*/

#ifndef waveAtmBoundaryLayerSuperposition_H
#define waveAtmBoundaryLayerSuperposition_H

#include "waveSuperposition.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class waveAtmBoundaryLayerSuperposition Declaration
\*---------------------------------------------------------------------------*/

class waveAtmBoundaryLayerSuperposition
:
    public waveSuperposition
{
    // Private member data

        //- The reference velocity
        const vector UGasRef_;

        //- The reference height
        const scalar hRef_;

        //- The minimum wave elevation
        const scalar hWaveMin_;

        //- The maximum wave elevation
        const scalar hWaveMax_;


public:

    //- Runtime type information
    TypeName("waveAtmBoundaryLayer");


    // Constructors

        //- Construct from a database
        waveAtmBoundaryLayerSuperposition(const objectRegistry& db);


    //- Destructor
    ~waveAtmBoundaryLayerSuperposition();


    // Member Functions

        //- Get the gas velocity at a given time and global positions
        virtual tmp<vectorField> UGas
        (
            const scalar t,
            const vectorField& p
        ) const;

        //- Inherit write from regIOobject
        using regIOobject::write;

        //- Write
        void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
