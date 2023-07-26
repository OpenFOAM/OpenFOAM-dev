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
    Foam::compressible::cavitationModels::Saito

Description
    Saito cavitation model.

    Reference:
    \verbatim
        Saito, Y., Takami, R., Nakamori, I., & Ikohagi, T. (2007).
        Numerical analysis of unsteady behavior of cloud cavitation
        around a NACA0015 foil. Computational Mechanics, 40, 85-96.
   \endverbatim

Usage:
    \table
        Property | Description                      | Required | Default value
        liquid   | Name of the liquid phase         | yes      |
        pSat     | Saturation vapor pressure        | yes      |
        Ca       | Interfacial area concentration coefficient [1/m] | yes      |
        Cv       | Vapourisation rate coefficient   | yes      |
        Cc       | Condensation rate coefficient    | yes      |
        alphaNuc | Nucleation site volume fraction  | yes      |
    \endtable

    Example:
    \verbatim
        model       Saito;

        liquid      liquid;

        pSat        2300;

        Ca          0.1;  // Interfacial area concentration coefficient [1/m]
        Cc          1;
        Cv          1;
        alphaNuc    0.001;
    \endverbatim

SourceFiles
    Saito.C

\*---------------------------------------------------------------------------*/

#ifndef Saito_H
#define Saito_H

#include "cavitationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace cavitationModels
{

/*---------------------------------------------------------------------------*\
                                 Class Saito
\*---------------------------------------------------------------------------*/

class Saito
:
    public cavitationModel
{
    // Private Data

        //- Interfacial area concentration coefficient [1/m]
        dimensionedScalar Ca_;

        //- Vapourisation rate coefficient
        dimensionedScalar Cv_;

        //- Condensation rate coefficient
        dimensionedScalar Cc_;

        //- Nucleation site volume fraction
        dimensionedScalar alphaNuc_;

        //- Zero with pressure dimensions
        dimensionedScalar p0_;


    // Private Member functions

        //- Return the function of temperature
        tmp<volScalarField::Internal> fT(const rhoFluidThermo& thermo) const;


public:

    //- Runtime type information
    TypeName("Saito");


    // Constructors

        //- Construct for phases
        Saito
        (
            const dictionary& dict,
            const compressibleTwoPhases& phases
        );


    //- Destructor
    virtual ~Saito()
    {}


    // Member Functions

        //- Return the mass condensation and vaporisation rates as a
        //  coefficient to multiply alphav for the condensation rate and a
        //  coefficient to multiply alphal for the vaporisation rate
        virtual Pair<tmp<volScalarField::Internal>> mDotcvAlphal() const;

        //- Return the mass condensation and vaporisation rates as coefficients
        //  to multiply (p - pSat)
        virtual Pair<tmp<volScalarField::Internal>> mDotcvP() const;

        //- Correct the Saito phaseChange model
        virtual void correct();

        //- Read the dictionary and update
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace cavitationModels
} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
