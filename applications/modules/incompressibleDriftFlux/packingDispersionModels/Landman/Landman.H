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
    Foam::packingDispersionModels::Landman

Description
    Landman packing dispersion model for industrial separations

Usage
    Example usage:
    \verbatim
        packingDispersionModel Landman;

        LandmanCoeffs
        {
            sigma0          2.11e-3;
            n               5;
            alphaGel        0.14;
        }
    \endverbatim

SourceFiles
    Landman.C

\*---------------------------------------------------------------------------*/

#ifndef Landman_H
#define Landman_H

#include "packingDispersionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace packingDispersionModels
{

/*---------------------------------------------------------------------------*\
                           Class Landman Declaration
\*---------------------------------------------------------------------------*/

class Landman
:
    public packingDispersionModel
{
    // Private Data

        dimensionedScalar sigma0_;
        dimensionedScalar n_;
        dimensionedScalar alphaGel_;


public:

    //- Runtime type information
    TypeName("Landman");


    // Constructors

        //- Construct from components
        Landman
        (
            const dictionary& dict,
            const relativeVelocityModel& relativeVelocity
        );


    //- Destructor
    ~Landman();


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
