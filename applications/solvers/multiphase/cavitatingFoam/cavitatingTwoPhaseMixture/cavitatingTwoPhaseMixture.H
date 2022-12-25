/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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
    Foam::cavitatingTwoPhaseMixture

Description
    Class to represent a mixture of two constant density phases

SourceFiles
    cavitatingTwoPhaseMixture.C

\*---------------------------------------------------------------------------*/

#ifndef cavitatingTwoPhaseMixture_H
#define cavitatingTwoPhaseMixture_H

#include "twoPhaseMixture.H"
#include "viscosityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class barotropicCompressibilityModel;

/*---------------------------------------------------------------------------*\
                Class cavitatingTwoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class cavitatingTwoPhaseMixture
:
    public twoPhaseMixture,
    public viscosity
{
    // Private Data

        //- Vapour phase fraction reference (alpha1 from twoPhaseMixture)
        volScalarField& alphav_;

        //- liquid phase fraction reference (alpha2 from twoPhaseMixture)
        volScalarField& alphal_;

        //- Vapour viscosity model
        autoPtr<viscosityModel> nuModelv_;

        //- Liquid viscosity model
        autoPtr<viscosityModel> nuModell_;

        //- Constant vapour density
        dimensionedScalar rhov_;

        //- Constant liquid density
        dimensionedScalar rhol_;

        //- Mixture viscosity
        volScalarField nu_;

        //- Thermophysical properties dictionary
        IOdictionary thermodynamicProperties_;

        //- Liquid compressibility
        dimensionedScalar psil_;

        //- Saturation liquid density
        dimensionedScalar rholSat_;

        //- Vapour compressibility
        dimensionedScalar psiv_;

        //- Saturation pressure
        dimensionedScalar pSat_;

        //- Saturation vapour density
        dimensionedScalar rhovSat_;

        //- rholSat_ - pSat_*psil_
        dimensionedScalar rhol0_;

        //- Minimum density for numerical stability
        dimensionedScalar rhoMin_;

        //- Static pressure field
        volScalarField p_;

        autoPtr<barotropicCompressibilityModel> psiModel_;

        //- Reference to the mixture compressibility field
        const volScalarField& psi_;

        //- Mixture density field
        volScalarField rho_;


public:

    TypeName("cavitatingTwoPhaseMixture");


    // Constructors

        //- Construct from a mesh
        cavitatingTwoPhaseMixture(const fvMesh& mesh);


    //- Destructor
    virtual ~cavitatingTwoPhaseMixture();


    // Member Functions

        //- Return the static pressure field
        volScalarField& p()
        {
            return p_;
        }

        //- Return the mixture density field
        virtual const volScalarField& rho() const
        {
            return rho_;
        }

        //- Return the mixture density field
        volScalarField& rho()
        {
            return rho_;
        }

        //- Return the mixture compressibility
        const volScalarField& psi() const
        {
            return psi_;
        }

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

        //- Correct the pressure from mixture density
        void correctPressure();

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
