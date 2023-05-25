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
    Foam::incompressibleTwoPhaseVoFMixture

Description
    Class to represent a mixture of two constant density phases

SourceFiles
    incompressibleTwoPhaseVoFMixture.C

\*---------------------------------------------------------------------------*/

#ifndef incompressibleTwoPhaseVoFMixture_H
#define incompressibleTwoPhaseVoFMixture_H

#include "twoPhaseVoFMixture.H"
#include "incompressibleTwoPhases.H"
#include "viscosityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class incompressibleTwoPhaseVoFMixture Declaration
\*---------------------------------------------------------------------------*/

class incompressibleTwoPhaseVoFMixture
:
    public twoPhaseVoFMixture,
    virtual public incompressibleTwoPhases,
    public viscosity
{
    // Private Data

        //- Viscosity model for phase 1
        autoPtr<viscosityModel> nuModel1_;

        //- Viscosity model for phase 2
        autoPtr<viscosityModel> nuModel2_;

        //- Constant density of phase 1
        dimensionedScalar rho1_;

        //- Constant density of phase 2
        dimensionedScalar rho2_;

        //- Mixture density
        volScalarField rho_;

        //- Mixture viscosity
        volScalarField nu_;


public:

    TypeName("incompressibleTwoPhaseVoFMixture");


    // Constructors

        //- Construct from a mesh
        incompressibleTwoPhaseVoFMixture(const fvMesh& mesh);


    //- Destructor
    virtual ~incompressibleTwoPhaseVoFMixture()
    {}


    // Member Functions

        //- Return const-access to phase1 viscosityModel
        const viscosityModel& nuModel1() const
        {
            return nuModel1_();
        }

        //- Return const-access to phase2 viscosityModel
        const viscosityModel& nuModel2() const
        {
            return nuModel2_();
        }

        //- Return const-access to phase1 density
        const dimensionedScalar& rho1() const
        {
            return rho1_;
        }

        //- Return const-access to phase2 density
        const dimensionedScalar& rho2() const
        {
            return rho2_;
        };

        //- Return the mixture density
        virtual const volScalarField& rho() const
        {
            return rho_;
        }

        //- Return the dynamic laminar viscosity
        tmp<volScalarField> mu() const;

        //- Return the face-interpolated dynamic laminar viscosity
        tmp<surfaceScalarField> muf() const;

        //- Return the kinematic laminar viscosity
        virtual tmp<volScalarField> nu() const
        {
            return nu_;
        }

        //- Return the laminar viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const
        {
            return nu_.boundaryField()[patchi];
        }

        //- Return the face-interpolated kinematic laminar viscosity
        tmp<surfaceScalarField> nuf() const;

        //- Correct the mixture density and laminar viscosity
        virtual void correct();

        //- Read base phaseProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
