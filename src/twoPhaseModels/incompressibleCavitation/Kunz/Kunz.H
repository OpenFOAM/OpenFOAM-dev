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
    Foam::cavitationModels::Kunz

Description
    Kunz cavitation model.

    Slightly modified so that the condensation term is switched off when the
    pressure is less than the saturation vapour pressure.  This change allows
    the condensation term to be formulated as a coefficient multiplying (p -
    p_sat) so that it can be included as an implicit term in the pressure
    equation.

    Reference:
    \verbatim
        Kunz, R. F., Boger, D. A., Stinebring, D. R., Chyczewski,
        T. S., Lindau, J. W., Gibeling, H. J., ... & Govindan, T. (2000).
        A preconditioned Navier–Stokes method for two-phase flows
        with application to cavitation prediction.
        Computers & Fluids, 29(8), 849-875.
    \endverbatim

Usage:
    \table
        Property | Description                     | Required | Default value
        liquid   | Name of the liquid phase             | yes      |
        pSat     | Saturation vapor pressure            | yes      |
        Uinf     | Free-stream velocity                 | yes      |
        tInf     | Vaporisation/condensation time scale | yes      |
        Cv       | Vapourisation rate coefficient       | yes      |
        Cc       | Condensation rate coefficient        | yes      |
    \endtable

    Example:
    \verbatim
        model       Kunz;

        liquid      liquid;

        pSat        2300;

        UInf        20;
        tInf        0.005; // L = 0.1 m
        Cv          1000;
        Cc          1000;
    \endverbatim

SourceFiles
    Kunz.C

\*---------------------------------------------------------------------------*/

#ifndef Kunz_H
#define Kunz_H

#include "cavitationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace cavitationModels
{

/*---------------------------------------------------------------------------*\
                                 Class Kunz
\*---------------------------------------------------------------------------*/

class Kunz
:
    public cavitationModel
{
    // Private Data

        //- Free-stream velocity
        dimensionedScalar UInf_;

        //- Vaporisation time scale
        dimensionedScalar tInf_;

        //- Vaporisation rate coefficient
        dimensionedScalar Cv_;

        //- Condensation rate coefficient
        dimensionedScalar Cc_;

        //- Zero with pressure dimensions
        dimensionedScalar p0_;

        dimensionedScalar mcCoeff_;

        dimensionedScalar mvCoeff_;


public:

    //- Runtime type information
    TypeName("Kunz");


    // Constructors

        //- Construct for phases
        Kunz
        (
            const dictionary& dict,
            const incompressibleTwoPhases& phases
        );


    //- Destructor
    virtual ~Kunz()
    {}


    // Member Functions

        //- Return the mass condensation and vaporisation rates as a
        //  coefficient to multiply alphav for the condensation rate and a
        //  coefficient to multiply alphal for the vaporisation rate
        virtual Pair<tmp<volScalarField::Internal>> mDotcvAlpha() const;

        //- Return the mass condensation and vaporisation rates as coefficients
        //  to multiply (p - pSat)
        virtual Pair<tmp<volScalarField::Internal>> mDotcvP() const;

        //- Correct the Kunz phaseChange model
        virtual void correct();

        //- Read the dictionary and update
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace cavitationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
