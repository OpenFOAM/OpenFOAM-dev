/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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
    Foam::relativeVelocityModels::MichaelsBolger

Description
    Michaels & Bolger relative velocity model

    Reference:
    \verbatim
        Michaels, A. S., & Bolger, J. C. (1962).
        Settling rates and sediment volumes
        of flocculated kaolin suspensions.
        Industrial & Engineering Chemistry Fundamentals, 1(1), 24-33.
    \endverbatim

Usage
    Example usage:
    \verbatim
        relativeVelocityModel MichaelsBolger;

        MichaelsBolgerCoeffs
        {
            a0          0;    // Extended Michaels & Bolger coefficient,
            a1          4.65; // Exponent
        }
    \endverbatim

SourceFiles
    MichaelsBolger.C

\*---------------------------------------------------------------------------*/

#ifndef MichaelsBolger_H
#define MichaelsBolger_H

#include "relativeVelocityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relativeVelocityModels
{

/*---------------------------------------------------------------------------*\
                           Class MichaelsBolger Declaration
\*---------------------------------------------------------------------------*/

class MichaelsBolger
:
    public relativeVelocityModel
{
    // Private Data

        //- Extended Michaels & Bolger coefficient, defaults to 0
        dimensionedScalar a0_;

        //- Exponent, defaults to 4.65
        dimensionedScalar a1_;

        //- Drift velocity coefficient
        dimensionedScalar Vc_;


public:

    //- Runtime type information
    TypeName("MichaelsBolger");


    // Constructors

        //- Construct from components
        MichaelsBolger
        (
            const dictionary& dict,
            const incompressibleDriftFluxMixture& mixture,
            const uniformDimensionedVectorField& g
        );


    //- Destructor
    ~MichaelsBolger();


    // Member Functions

        //- Update the diffusion velocity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace relativeVelocityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
