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
    Foam::packingDispersionModels::Green

Description
    Green packing dispersion model for industrial separations

    Reference:
    \verbatim
        Green, M. D., Eberl, M., & Landman, K. A. (1996).
        Compressive yield stress of flocculated suspensions:
        determination via experiment.
        AIChE journal, 42(8), 2308-2318.
   \endverbatim

Usage
    Example usage:
    \verbatim
        packingDispersionModel Green;

        GreenCoeffs
        {
            sigma0          2.11e-3;
            n               9;
            alphaGel        0.1;
        }
    \endverbatim

SourceFiles
    Green.C

\*---------------------------------------------------------------------------*/

#ifndef Green_H
#define Green_H

#include "packingDispersionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace packingDispersionModels
{

/*---------------------------------------------------------------------------*\
                           Class Green Declaration
\*---------------------------------------------------------------------------*/

class Green
:
    public packingDispersionModel
{
    // Private Data

        dimensionedScalar sigma0_;
        dimensionedScalar n_;
        dimensionedScalar alphaGel_;


public:

    //- Runtime type information
    TypeName("Green");


    // Constructors

        //- Construct from components
        Green
        (
            const dictionary& dict,
            const relativeVelocityModel& relativeVelocity
        );


    //- Destructor
    ~Green();


    // Member Functions

        //- Return the derivative of the packing stress w.r.t. phase-fraction
        virtual tmp<volScalarField> sigmaPrime() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace packingDispersionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
