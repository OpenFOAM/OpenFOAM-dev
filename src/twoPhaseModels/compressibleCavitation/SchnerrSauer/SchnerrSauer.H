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
    Foam::compressible::cavitationModels::SchnerrSauer

Description
    SchnerrSauer cavitation model.

    Reference:
    \verbatim
        Schnerr, G. H., & Sauer, J. (2001, May).
        Physical and numerical modeling of unsteady cavitation dynamics.
        In Fourth international conference on multiphase flow (Vol. 1).
        New Orleans, LO, USA: ICMF New Orleans.
    \endverbatim

Usage:
    \table
        Property | Description                    | Required | Default value
        liquid   | Name of the liquid phase       | yes      |
        pSat     | Saturation vapor pressure      | yes      |
        Uinf     | Free-stream velocity           | yes      |
        n        | Bubble number density          | yes      |
        dNuc     | Nucleation site diameter       | yes      |
        Cv       | Vapourisation rate coefficient | yes      |
        Cc       | Condensation rate coefficient  | yes      |
    \endtable

    Example:
    \verbatim
        model       SchnerrSauer;

        liquid      liquid;

        pSat        2300;

        n           1.6e+13;
        dNuc        2e-6;
        Cv          1;
        Cc          1;
    \endverbatim

SourceFiles
    SchnerrSauer.C

\*---------------------------------------------------------------------------*/

#ifndef SchnerrSauer_H
#define SchnerrSauer_H

#include "cavitationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace cavitationModels
{

/*---------------------------------------------------------------------------*\
                              Class SchnerrSauer
\*---------------------------------------------------------------------------*/

class SchnerrSauer
:
    public cavitationModel
{
    // Private Data

        //- Bubble number density
        dimensionedScalar n_;

        //- Nucleation site diameter
        dimensionedScalar dNuc_;

        //- Vaporisation rate coefficient
        dimensionedScalar Cv_;

        //- Condensation rate coefficient
        dimensionedScalar Cc_;

        //- Zero with pressure dimensions
        dimensionedScalar p0_;


    // Private Member functions

        //- Nucleation site volume-fraction
        dimensionedScalar alphaNuc() const;

        //- Reciprocal bubble radius
        tmp<volScalarField::Internal> rRb
        (
            const volScalarField::Internal& limitedAlphal
        ) const;

        //- Part of the condensation and vaporisation rates
        tmp<volScalarField::Internal> pCoeff
        (
            const volScalarField::Internal& p,
            const volScalarField::Internal& pSat
        ) const;


public:

    //- Runtime type information
    TypeName("SchnerrSauer");


    // Constructors

        //- Construct for phases
        SchnerrSauer
        (
            const dictionary& dict,
            const compressibleTwoPhases& phases
        );


    //- Destructor
    virtual ~SchnerrSauer()
    {}


    // Member Functions

        //- Return the mass condensation and vaporisation rates as a
        //  coefficient to multiply alphav for the condensation rate and a
        //  coefficient to multiply alphal for the vaporisation rate
        virtual Pair<tmp<volScalarField::Internal>> mDotcvAlphal() const;

        //- Return the mass condensation and vaporisation rates as coefficients
        //  to multiply (p - pSat)
        virtual Pair<tmp<volScalarField::Internal>> mDotcvP() const;

        //- Correct the SchnerrSauer phaseChange model
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
