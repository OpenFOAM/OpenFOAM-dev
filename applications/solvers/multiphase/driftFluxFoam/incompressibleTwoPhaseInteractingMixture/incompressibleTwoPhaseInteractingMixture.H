/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2022 OpenFOAM Foundation
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
    Foam::incompressibleTwoPhaseInteractingMixture

Description
    A two-phase incompressible transportModel for interacting phases
    requiring the direct evaluation of the mixture viscosity,
    e.g. activated sludge or slurry.

SourceFiles
    incompressibleTwoPhaseInteractingMixture.C

\*---------------------------------------------------------------------------*/

#ifndef incompressibleTwoPhaseInteractingMixture_H
#define incompressibleTwoPhaseInteractingMixture_H

#include "twoPhaseMixture.H"
#include "viscosityModel.H"
#include "uniformDimensionedFields.H"
#include "IOMRFZoneList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mixtureViscosityModel;
class relativeVelocityModel;

/*---------------------------------------------------------------------------*\
          Class incompressibleTwoPhaseInteractingMixture Declaration
\*---------------------------------------------------------------------------*/

class incompressibleTwoPhaseInteractingMixture
:
    public twoPhaseMixture,
    public viscosity
{
    // Private data

        volVectorField& U_;

        autoPtr<mixtureViscosityModel> muModel_;
        autoPtr<viscosityModel> nucModel_;

        dimensionedScalar rhod_;
        dimensionedScalar rhoc_;

        //- Optional diameter of the dispersed phase particles
        dimensionedScalar dd_;

        //- Optional maximum dispersed phase-fraction (e.g. packing limit)
        scalar alphaMax_;

        volScalarField rho_;

        volScalarField mu_;

        //- Optional MRF zones
        IOMRFZoneList MRF_;

        //- Dispersed phase relative velocity model
        autoPtr<relativeVelocityModel> UdmModel_;


public:

    TypeName("incompressibleTwoPhaseInteractingMixture");


    // Constructors

        //- Construct from components
        incompressibleTwoPhaseInteractingMixture
        (
            volVectorField& U,
            const surfaceScalarField& phi,
            const uniformDimensionedVectorField& g
        );


    //- Destructor
    virtual ~incompressibleTwoPhaseInteractingMixture();


    // Member Functions

        //- Return const-access to the continuous phase-fraction
        const volScalarField& alphac() const;

        //- Return const-access to the dispersed phase-fraction
        const volScalarField& alphad() const;

        //- Return const-access to the mixture viscosityModel
        const mixtureViscosityModel& muModel() const;

        //- Return const-access to the continuous-phase viscosityModel
        const viscosityModel& nucModel() const;

        //- Return const-access to the dispersed-phase density
        const dimensionedScalar& rhod() const;

        //- Return const-access to continuous-phase density
        const dimensionedScalar& rhoc() const;

        //- Return the diameter of the dispersed-phase particles
        const dimensionedScalar& dd() const;

        //- Optional maximum phase-fraction (e.g. packing limit)
        //  Defaults to 1
        scalar alphaMax() const;

        //- Return const-access to the mixture velocity
        const volVectorField& U() const;

        //- Return the mixture density
        virtual const volScalarField& rho() const;

        //- Return the dynamic mixture viscosity
        tmp<volScalarField> mu() const;

        //- Return the dynamic mixture viscosity for patch
        virtual tmp<scalarField> mu(const label patchi) const;

        //- Return the mixture viscosity
        virtual tmp<volScalarField> nu() const;

        //- Return the mixture viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const;

        //- Return MRF zones
        const IOMRFZoneList& MRF() const;

        //- Return the diffusion velocity of the dispersed phase
        const volVectorField& Udm() const;

        //- Return the div stress tensor due to the phase transport
        tmp<volVectorField> divTauDm() const;

        //- Correct the laminar viscosity
        virtual void correct();

        //- Read base phaseProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
