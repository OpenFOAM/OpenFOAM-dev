/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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
    Foam::PhaseCompressibleTurbulenceModel

Description
    Templated abstract base class for multiphase compressible
    turbulence models.

SourceFiles
    PhaseCompressibleTurbulenceModel.C

\*---------------------------------------------------------------------------*/

#ifndef PhaseCompressibleTurbulenceModel_H
#define PhaseCompressibleTurbulenceModel_H

#include "TurbulenceModel.H"
#include "compressibleTurbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class PhaseCompressibleTurbulenceModel Declaration
\*---------------------------------------------------------------------------*/

template<class TransportModel>
class PhaseCompressibleTurbulenceModel
:
    public TurbulenceModel
    <
        volScalarField,
        volScalarField,
        compressibleTurbulenceModel,
        TransportModel
    >
{

public:

    typedef volScalarField alphaField;
    typedef volScalarField rhoField;
    typedef TransportModel transportModel;


    // Constructors

        //- Construct
        PhaseCompressibleTurbulenceModel
        (
            const word& type,
            const alphaField& alpha,
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& trasport,
            const word& propertiesName
        );


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<PhaseCompressibleTurbulenceModel> New
        (
            const alphaField& alpha,
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& trasportModel,
            const word& propertiesName = turbulenceModel::propertiesName
        );


    //- Destructor
    virtual ~PhaseCompressibleTurbulenceModel()
    {}


    // Member Functions

        //- Return the laminar dynamic viscosity
        virtual tmp<volScalarField> mu() const
        {
            return this->transport_.mu();
        }

        //- Return the laminar dynamic viscosity on patch
        virtual tmp<scalarField> mu(const label patchi) const
        {
            return this->transport_.mu(patchi);
        }

        //- Return the turbulence dynamic viscosity
        virtual tmp<volScalarField> mut() const
        {
            return this->rho_*this->nut();
        }

        //- Return the turbulence dynamic viscosity on patch
        virtual tmp<scalarField> mut(const label patchi) const
        {
            return this->rho_.boundaryField()[patchi]*this->nut(patchi);
        }

        //- Return the effective dynamic viscosity
        virtual tmp<volScalarField> muEff() const
        {
            return mut() + mu();
        }

        //- Return the effective dynamic viscosity on patch
        virtual tmp<scalarField> muEff(const label patchi) const
        {
            return mut(patchi) + mu(patchi);
        }

        //- Return the phase-pressure'
        // (derivative of phase-pressure w.r.t. phase-fraction)
        virtual tmp<volScalarField> pPrime() const;

        //- Return the face-phase-pressure'
        // (derivative of phase-pressure w.r.t. phase-fraction)
        virtual tmp<surfaceScalarField> pPrimef() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PhaseCompressibleTurbulenceModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
