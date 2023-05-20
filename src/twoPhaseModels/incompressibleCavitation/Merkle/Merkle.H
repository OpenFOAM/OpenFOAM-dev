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
    Foam::cavitationModels::Merkle

Description
    Merkle cavitation model.

    Reference:
    \verbatim
        Markle, C. L. (1998).
        Computational modeling of the dynamics of sheet cavitation.
        In Proc. 3rd International Symposium on Cavitation
        (Vol. 2, pp. 307-311).
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
        model       Merkle;

        liquid      liquid;

        pSat        2300;

        UInf        20;
        tInf        0.005;  // L = 0.1 m
        Cv          1e-3;
        Cc          80;
    \endverbatim

SourceFiles
    Merkle.C

\*---------------------------------------------------------------------------*/

#ifndef Merkle_H
#define Merkle_H

#include "cavitationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace cavitationModels
{

/*---------------------------------------------------------------------------*\
                                 Class Merkle
\*---------------------------------------------------------------------------*/

class Merkle
:
    public cavitationModel
{
    // Private Data

        //- Free-stream velocity
        dimensionedScalar UInf_;

        //- Vaporisation/condensation time scale
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
    TypeName("Merkle");


    // Constructors

        //- Construct for phases
        Merkle
        (
            const dictionary& dict,
            const incompressibleTwoPhases& phases
        );


    //- Destructor
    virtual ~Merkle()
    {}


    // Member Functions

        //- Return the mass condensation and vaporisation rates as a
        //  coefficient to multiply alphav for the condensation rate and a
        //  coefficient to multiply alphal for the vaporisation rate
        virtual Pair<tmp<volScalarField::Internal>> mDotcvAlpha() const;

        //- Return the mass condensation and vaporisation rates as coefficients
        //  to multiply (p - pSat)
        virtual Pair<tmp<volScalarField::Internal>> mDotcvP() const;

        //- Correct the Merkle phaseChange model
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
