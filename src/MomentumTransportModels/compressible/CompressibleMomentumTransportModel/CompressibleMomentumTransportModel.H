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
    Foam::CompressibleMomentumTransportModel

Description
    Templated abstract base class for single-phase compressible
    turbulence models.

SourceFiles
    CompressibleMomentumTransportModel.C

\*---------------------------------------------------------------------------*/

#ifndef CompressibleMomentumTransportModel_H
#define CompressibleMomentumTransportModel_H

#include "MomentumTransportModel.H"
#include "compressibleMomentumTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class CompressibleMomentumTransportModel Declaration
\*---------------------------------------------------------------------------*/

template<class TransportModel>
class CompressibleMomentumTransportModel
:
    public MomentumTransportModel
    <
        geometricOneField,
        volScalarField,
        compressibleMomentumTransportModel,
        TransportModel
    >
{

public:

    typedef geometricOneField alphaField;
    typedef volScalarField rhoField;
    typedef TransportModel transportModel;


    // Constructors

        //- Construct
        CompressibleMomentumTransportModel
        (
            const word& type,
            const geometricOneField& alpha,
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport
        );


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<CompressibleMomentumTransportModel> New
        (
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const transportModel& transportModel
        );


    //- Destructor
    virtual ~CompressibleMomentumTransportModel()
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

        //- Return the laminar viscosity
        virtual tmp<volScalarField> nu() const
        {
            return this->transport_.mu()/this->rho_;
        }

        //- Return the laminar viscosity on patchi
        virtual tmp<scalarField> nu(const label patchi) const
        {
            return
                this->transport_.mu(patchi)
               /this->rho_.boundaryField()[patchi];
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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CompressibleMomentumTransportModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
