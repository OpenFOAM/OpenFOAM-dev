/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2021 OpenFOAM Foundation
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

#include "viscosityModel.H"
#include "mixtureViscosityModel.H"
#include "twoPhaseMixture.H"
#include "viscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class incompressibleTwoPhaseInteractingMixture Declaration
\*---------------------------------------------------------------------------*/

class incompressibleTwoPhaseInteractingMixture
:
    public twoPhaseMixture,
    public viscosity
{
    // Private data

        autoPtr<mixtureViscosityModel> muModel_;
        autoPtr<viscosityModel> nucModel_;

        dimensionedScalar rhod_;
        dimensionedScalar rhoc_;

        //- Optional diameter of the dispersed phase particles
        dimensionedScalar dd_;

        //- Optional maximum dispersed phase-fraction (e.g. packing limit)
        scalar alphaMax_;

        const volVectorField& U_;

        volScalarField mu_;


public:

    TypeName("incompressibleTwoPhaseInteractingMixture");


    // Constructors

        //- Construct from components
        incompressibleTwoPhaseInteractingMixture
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~incompressibleTwoPhaseInteractingMixture()
    {}


    // Member Functions

        //- Return const-access to the mixture viscosityModel
        const mixtureViscosityModel& muModel() const
        {
            return muModel_();
        }

        //- Return const-access to the continuous-phase viscosityModel
        const viscosityModel& nucModel() const
        {
            return nucModel_();
        }

        //- Return const-access to the dispersed-phase density
        const dimensionedScalar& rhod() const
        {
            return rhod_;
        }

        //- Return const-access to continuous-phase density
        const dimensionedScalar& rhoc() const
        {
            return rhoc_;
        };

        //- Return the diameter of the dispersed-phase particles
        const dimensionedScalar& dd() const
        {
            return dd_;
        }

        //- Optional maximum phase-fraction (e.g. packing limit)
        //  Defaults to 1
        scalar alphaMax() const
        {
            return alphaMax_;
        }

        //- Return const-access to the mixture velocity
        const volVectorField& U() const
        {
            return U_;
        }

        //- Return the dynamic mixture viscosity
        tmp<volScalarField> mu() const
        {
            return mu_;
        }

        //- Return the dynamic mixture viscosity for patch
        virtual tmp<scalarField> mu(const label patchi) const
        {
            return mu_.boundaryField()[patchi];
        }

        //- Return the mixture density
        virtual tmp<volScalarField> rho() const
        {
            return alpha1()*rhod_ + alpha2()*rhoc_;
        }

        //- Return the mixture density for patch
        virtual tmp<scalarField> rho(const label patchi) const
        {
            return
                alpha1().boundaryField()[patchi]*rhod_.value()
              + alpha2().boundaryField()[patchi]*rhoc_.value();
        }

        //- Return the mixture viscosity
        virtual tmp<volScalarField> nu() const
        {
            return mu_/rho();
        }

        //- Return the mixture viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const
        {
            return mu_.boundaryField()[patchi]/rho(patchi);
        }

        //- Correct the laminar viscosity
        virtual void correct()
        {
            mu_ = muModel_->mu(rhoc_*nucModel_->nu(), U_);
        }

        //- Read base phaseProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
