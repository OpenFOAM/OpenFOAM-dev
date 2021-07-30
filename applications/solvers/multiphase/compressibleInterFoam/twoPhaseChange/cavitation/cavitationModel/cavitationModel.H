/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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
    Foam::twoPhaseChangeModels::cavitationModel

Description
    Abstract base class for cavitation models

SourceFiles
    cavitationModel.C

\*---------------------------------------------------------------------------*/

#ifndef cavitationModel_H
#define cavitationModel_H

#include "twoPhaseChangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace twoPhaseChangeModels
{

/*---------------------------------------------------------------------------*\
                                 Class cavitationModel
\*---------------------------------------------------------------------------*/

class cavitationModel
:
    public twoPhaseChangeModel
{
    // Private data

        //- Saturation vapour pressure
        dimensionedScalar pSat_;


public:

    //- Runtime type information
    TypeName("cavitation");


    // Constructors

        //- Construct for mixture
        cavitationModel
        (
            const word& type,
            const compressibleTwoPhaseMixture& mixture
        );


    //- Destructor
    virtual ~cavitationModel()
    {}


    // Member Functions

        //- Return const-access to the saturation vapour pressure
        const dimensionedScalar& pSat() const
        {
            return pSat_;
        }

        //- Return the mass condensation and vaporisation rates as a
        //  coefficient to multiply (1 - alphal) for the condensation rate
        //  and a coefficient to multiply  alphal for the vaporisation rate
        virtual Pair<tmp<volScalarField::Internal>> mDotAlphal() const = 0;

        //- Return the mass condensation and vaporisation rates as coefficients
        //  to multiply (p - pSat)
        virtual Pair<tmp<volScalarField::Internal>> mDotP() const = 0;

        //- Return the cavitation explicit and implicit sources
        //  for the phase-fraction equation
        virtual Pair<tmp<volScalarField::Internal>> Salpha
        (
            volScalarField& alpha
        ) const;

        //- Return the cavitation source matrix
        //  for the p_rgh pressure equation
        virtual tmp<fvScalarMatrix> Sp_rgh
        (
            const volScalarField& rho,
            const volScalarField& gh,
            volScalarField& p_rgh
        ) const;

        //- Return the cavitation source matrix for the momentum equation
        virtual tmp<fvVectorMatrix> SU
        (
            const volScalarField& rho,
            const surfaceScalarField& rhoPhi,
            volVectorField& U
        ) const;

        //- Read the phaseProperties dictionary and update
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace twoPhaseChangeModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
