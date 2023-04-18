/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2023 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mixtureViscosityModel;

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

        //- Continuous phase viscosity model
        autoPtr<viscosityModel> nucModel_;

        //- Mixture viscosity model
        autoPtr<mixtureViscosityModel> muModel_;

        //- Continuous phase density
        dimensionedScalar rhoc_;

        //- Dispersed phase density
        dimensionedScalar rhod_;

        //- Optional maximum dispersed phase-fraction (e.g. packing limit)
        scalar alphaMax_;

        //- Mixture density
        volScalarField rho_;

        //- Mixture kinematic viscosity
        volScalarField nu_;


public:

    TypeName("incompressibleTwoPhaseInteractingMixture");


    // Constructors

        //- Construct from components
        incompressibleTwoPhaseInteractingMixture
        (
            volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~incompressibleTwoPhaseInteractingMixture();


    // Member Functions

        //- Return const-access to the continuous phase-fraction
        const volScalarField& alphac() const;

        //- Return const-access to the dispersed phase-fraction
        const volScalarField& alphad() const;

        //- Return const-access to continuous-phase density
        const dimensionedScalar& rhoc() const;

        //- Return const-access to the dispersed-phase density
        const dimensionedScalar& rhod() const;

        //- Optional maximum phase-fraction (e.g. packing limit)
        //  Defaults to 1
        scalar alphaMax() const;

        //- Return const-access to the mixture velocity
        const volVectorField& U() const;

        //- Return the mixture density
        virtual const volScalarField& rho() const;

        //- Return the mixture viscosity
        virtual tmp<volScalarField> nu() const;

        //- Return the mixture viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const;

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
