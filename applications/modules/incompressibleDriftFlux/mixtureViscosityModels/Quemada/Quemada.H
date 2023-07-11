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
    Foam::mixtureViscosityModels::Quemada

Description
     Quemada viscosity model for colloidal dispersions.

     References:
     \verbatim
         Quemada, D. (1998).
         Rheological modelling of complex fluids. I.
         The concept of effective volume fraction revisited.
         The European Physical Journal-Applied Physics, 1(1), 119-127.
    \endverbatim

Usage
    Example usage:
    \verbatim
        viscosityModel  Quemada;

        rho         1996; // Dispersed phase density

        q           2;    // Exponent, defaults to 2
        muMax       1e-2; // Maximum viscosity (for numerical stability)

    \endverbatim

SourceFiles
    Quemada.C

\*---------------------------------------------------------------------------*/

#ifndef Quemada_H
#define Quemada_H

#include "mixtureViscosityModel.H"
#include "dimensionedScalar.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace mixtureViscosityModels
{

/*---------------------------------------------------------------------------*\
                           Class Quemada Declaration
\*---------------------------------------------------------------------------*/

class Quemada
:
    public mixtureViscosityModel
{
    // Private data

        //- Exponent (defaults to 2)
        scalar q_;

        //- Maximum viscosity
        dimensionedScalar muMax_;


public:

    //- Runtime type information
    TypeName("Quemada");


    // Constructors

        //- Construct from mixture
        Quemada(const incompressibleDriftFluxMixture& mixture);


    //- Destructor
    virtual ~Quemada()
    {}


    // Member Functions

        //- Return the mixture viscosity
        //  given the viscosity of the continuous phase
        virtual tmp<volScalarField> mu
        (
            const volScalarField& muc,
            const volVectorField& U
        ) const;

        //- Read phaseProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace mixtureViscosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
