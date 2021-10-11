/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2021 OpenFOAM Foundation
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
    Foam::compressibleTwoPhaseMixture

Description
    Thermophysical properties class for a mixture of two rhoThermo including
    interfaceProperties for VoF for use in compressibleInterFoam.

SourceFiles
    compressibleTwoPhaseMixtureI.H
    compressibleTwoPhaseMixture.C
    compressibleTwoPhaseMixtureIO.C

\*---------------------------------------------------------------------------*/

#ifndef compressibleTwoPhaseMixture_H
#define compressibleTwoPhaseMixture_H

#include "rhoThermo.H"
#include "twoPhaseMixture.H"
#include "interfaceProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class compressibleTwoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class compressibleTwoPhaseMixture
:
    public twoPhaseMixture,
    public interfaceProperties,
    public viscosity
{
    // Private Data

        //- Switch to choose between solving for internal energy
        //  or total internal energy which is the default
        Switch totalInternalEnergy_;

        //- Pressure
        volScalarField p_;

        //- Mixture temperature
        volScalarField T_;

        //- Thermo-package of phase 1
        autoPtr<rhoThermo> thermo1_;

        //- Thermo-package of phase 2
        autoPtr<rhoThermo> thermo2_;

        //- Mixture density
        volScalarField rho_;

        //- Mass-fraction of phase 1
        volScalarField Alpha1_;

        //- Mass-fraction of phase 2
        volScalarField Alpha2_;


public:

    //- Runtime type information
    TypeName("compressibleTwoPhaseMixture");


    // Constructors

        //- Construct from components
        compressibleTwoPhaseMixture
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~compressibleTwoPhaseMixture();


    // Member Functions

        //- Return true to solve for total internal energy
        //  return false to solve for internal energy
        bool totalInternalEnergy() const
        {
            return totalInternalEnergy_;
        }

        //- Return pressure [Pa]
        volScalarField& p()
        {
            return p_;
        }

        //- Return mixture temperature [K]
        volScalarField& T()
        {
            return T_;
        }

        //- Return the thermo for phase 1
        const rhoThermo& thermo1() const
        {
            return thermo1_();
        }

        //- Return the thermo for phase 2
        const rhoThermo& thermo2() const
        {
            return thermo2_();
        }

        //- Return the thermo for phase 1
        rhoThermo& thermo1()
        {
            return thermo1_();
        }

        //- Return the thermo for phase 2
        rhoThermo& thermo2()
        {
            return thermo2_();
        }

        //- Return the density of phase 1
        const volScalarField& rho1() const
        {
            return thermo1_->rho();
        }

        //- Return the density of phase 2
        const volScalarField& rho2() const
        {
            return thermo2_->rho();
        }

        //- Return mixture density [kg/m^3]
        const volScalarField& rho() const
        {
            return rho_;
        }

        //- Correct the thermodynamics of each phase
        virtual void correctThermo();

        //- Update mixture properties
        virtual void correct();

        //- Kinematic viscosity of mixture [m^2/s]
        virtual tmp<volScalarField> nu() const;

        //- Kinematic viscosity of mixture for patch [m^2/s]
        virtual tmp<scalarField> nu(const label patchi) const;


    // IO

        //- Read base phaseProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
