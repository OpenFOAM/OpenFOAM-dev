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
    Foam::incompressibleDriftFluxMixture

Description
    Class to represent a mixture of two constant density phases

SourceFiles
    incompressibleDriftFluxMixture.C

\*---------------------------------------------------------------------------*/

#ifndef incompressibleDriftFluxMixture_H
#define incompressibleDriftFluxMixture_H

#include "twoPhaseVoFMixture.H"
#include "incompressibleTwoPhases.H"
#include "viscosityModel.H"
#include "mixtureViscosityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class incompressibleDriftFluxMixture Declaration
\*---------------------------------------------------------------------------*/

class incompressibleDriftFluxMixture
:
    public twoPhaseVoFMixture,
    virtual public incompressibleTwoPhases,
    public viscosity
{
    // Private Data

        //- Continuous phase viscosity model
        autoPtr<viscosityModel> nucModel_;

        //- Mixture viscosity model
        autoPtr<mixtureViscosityModel> muModel_;

        //- Dispersed phase density
        dimensionedScalar rhod_;

        //- Continuous phase density
        dimensionedScalar rhoc_;

        //- Optional maximum dispersed phase-fraction (e.g. packing limit)
        scalar alphaMax_;

        //- Mixture density
        volScalarField rho_;

        //- Mixture kinematic viscosity
        volScalarField nu_;

        //- Pointer to the mixture velocity field
        //  Set by the initialise(U) function
        const volVectorField* Uptr_;


public:

    TypeName("incompressibleDriftFluxMixture");


    // Constructors

        //- Construct from a mesh
        incompressibleDriftFluxMixture(const fvMesh& mesh);


    //- Destructor
    virtual ~incompressibleDriftFluxMixture()
    {}


    // Member Functions

        //- Return const-access to the dispersed phase-fraction
        const volScalarField& alphad() const
        {
            return alpha1();
        }

        //- Return const-access to the continuous phase-fraction
        const volScalarField& alphac() const
        {
            return alpha2();
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
        }

        //- Return const-access to the dispersed-phase density
        const dimensionedScalar& rho1() const
        {
            return rhod_;
        }

        //- Return const-access to continuous-phase density
        const dimensionedScalar& rho2() const
        {
            return rhoc_;
        }

        //- Optional maximum phase-fraction (e.g. packing limit)
        //  Defaults to 1
        scalar alphaMax() const
        {
            return alphaMax_;
        }

        //- Return the mixture density
        virtual const volScalarField& rho() const
        {
            return rho_;
        }

        //- Return the mixture viscosity
        virtual tmp<volScalarField> nu() const
        {
            return nu_;
        }

        //- Return the mixture viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const
        {
            return nu_.boundaryField()[patchi];
        }

        //- Return the mixture velocity field
        const volVectorField& U() const
        {
            return *Uptr_;
        }

        //- Correct the laminar viscosity
        incompressibleDriftFluxMixture& initialise(const volVectorField& U);

        //- Correct the density and mixture viscosity
        virtual void correct();

        //- Read base phaseProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
