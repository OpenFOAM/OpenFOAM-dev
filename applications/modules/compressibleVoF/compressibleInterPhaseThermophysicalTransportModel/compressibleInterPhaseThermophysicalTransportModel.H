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
    Foam::compressibleInterPhaseThermophysicalTransportModel

Description
    Compressible two-phase VoF thermophysicalTransportModel

    Provides the mixture thermal conductivity and diffusivity for
    Fourier's temperature gradient heat flux model.

SourceFiles
    compressibleInterPhaseThermophysicalTransportModel.C

\*---------------------------------------------------------------------------*/

#ifndef compressibleInterPhaseThermophysicalTransportModel_H
#define compressibleInterPhaseThermophysicalTransportModel_H

#include "thermophysicalTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class compressibleInterPhaseTransportModel;

/*---------------------------------------------------------------------------*\
            Class compressibleInterPhaseThermophysicalTransportModel Declaration
\*---------------------------------------------------------------------------*/

class compressibleInterPhaseThermophysicalTransportModel
:
    public thermophysicalTransportModel
{
    // Private Data

        //- compressibleInterPhaseTransportModel
        const compressibleInterPhaseTransportModel& momentumTransport_;


public:

    // Constructors

        //- Construct from components
        compressibleInterPhaseThermophysicalTransportModel
        (
            const compressibleInterPhaseTransportModel& momentumTransport
        );

        //- Disallow default bitwise copy construction
        compressibleInterPhaseThermophysicalTransportModel
        (
            const compressibleInterPhaseThermophysicalTransportModel&
        ) = delete;


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read();

        //- Const access to the coefficients dictionary
        virtual const dictionary& coeffDict() const;

        //- Effective thermal turbulent conductivity
        //  of mixture [W/m/K]
        virtual tmp<volScalarField> kappaEff() const ;

        //- Effective thermal turbulent conductivity
        //  of mixture for patch [W/m/K]
        virtual tmp<scalarField> kappaEff(const label patchi) const;

        //- Return the effective temperature transport coefficient
        //  derived from the phase internal energy equations i.e. <kappa/Cv>
        tmp<volScalarField> alphaEff() const;

        //- Return the heat flux [W/m^2]
        virtual tmp<surfaceScalarField> q() const;

        //- Return the patch heat flux correction [W/m^2]
        //  For isotropic or patch-aligned thermal conductivity qCorr is null
        virtual tmp<scalarField> qCorr(const label patchi) const;

        //- Return the source term for the energy equation
        virtual tmp<fvScalarMatrix> divq(volScalarField& he) const;

        //- Predict the transport coefficients if possible
        //  without solving thermophysical transport model equations
        virtual void predict();

        //- Solve the thermophysical transport model equations
        //  and correct the transport coefficients
        virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=
        (
            const compressibleInterPhaseThermophysicalTransportModel&
        ) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
