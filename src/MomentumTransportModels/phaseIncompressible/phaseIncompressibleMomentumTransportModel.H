/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2022 OpenFOAM Foundation
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
    Foam::phaseIncompressibleMomentumTransportModel

Description
    Templated abstract base class for multiphase incompressible
    turbulence models.

SourceFiles
    phaseIncompressibleMomentumTransportModel.C

\*---------------------------------------------------------------------------*/

#ifndef phaseIncompressibleMomentumTransportModel_H
#define phaseIncompressibleMomentumTransportModel_H

#include "incompressibleMomentumTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class phaseIncompressibleMomentumTransportModel Declaration
\*---------------------------------------------------------------------------*/

class phaseIncompressibleMomentumTransportModel
:
    public incompressibleMomentumTransportModel
{

public:

    typedef volScalarField alphaField;
    typedef geometricOneField rhoField;


protected:

    // Protected data

        const alphaField& alpha_;


public:

    // Declare run-time constructor selection table

        declareRunTimeNewSelectionTable
        (
            autoPtr,
            phaseIncompressibleMomentumTransportModel,
            dictionary,
            (
                const alphaField& alpha,
                const rhoField& rho,
                const volVectorField& U,
                const surfaceScalarField& alphaRhoPhi,
                const surfaceScalarField& phi,
                const viscosity& viscosity
            ),
            (alpha, rho, U, alphaRhoPhi, phi, viscosity)
        );


    // Constructors

        //- Construct
        phaseIncompressibleMomentumTransportModel
        (
            const word& type,
            const alphaField& alpha,
            const geometricOneField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const viscosity& viscosity
        );


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<phaseIncompressibleMomentumTransportModel> New
        (
            const alphaField& alpha,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const viscosity& viscosity
        );


    //- Destructor
    virtual ~phaseIncompressibleMomentumTransportModel()
    {}


    // Member Functions

        //- Access function to phase fraction
        const alphaField& alpha() const
        {
            return alpha_;
        }

        //- Return the phase-pressure'
        // (derivative of phase-pressure w.r.t. phase-fraction)
        virtual tmp<volScalarField> pPrime() const;

        //- Return the face-phase-pressure'
        // (derivative of phase-pressure w.r.t. phase-fraction)
        virtual tmp<surfaceScalarField> pPrimef() const;

        //- Return the effective stress tensor
        virtual tmp<volSymmTensorField> devSigma() const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevSigma(volVectorField& U) const;

        //- Return the effective stress tensor
        virtual tmp<volSymmTensorField> devTau() const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevTau(volVectorField& U) const;
};


namespace phaseIncompressible
{
    typedef phaseIncompressibleMomentumTransportModel momentumTransportModel;

    template<class PhaseIncompressibleMomentumTransportModel>
    autoPtr<PhaseIncompressibleMomentumTransportModel> New
    (
        const volScalarField& alpha,
        const volVectorField& U,
        const surfaceScalarField& alphaRhoPhi,
        const surfaceScalarField& phi,
        const viscosity& viscosity
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "phaseIncompressibleMomentumTransportModelTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
