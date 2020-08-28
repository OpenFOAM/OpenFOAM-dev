/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2020 OpenFOAM Foundation
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
    Foam::IncompressibleMomentumTransportModel

Description
    Templated abstract base class for single-phase incompressible
    turbulence models.

SourceFiles
    IncompressibleMomentumTransportModel.C

\*---------------------------------------------------------------------------*/

#ifndef IncompressibleMomentumTransportModel_H
#define IncompressibleMomentumTransportModel_H

#include "MomentumTransportModel.H"
#include "incompressibleMomentumTransportModel.H"
#include "fvMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class IncompressibleMomentumTransportModel Declaration
\*---------------------------------------------------------------------------*/

template<class TransportModel>
class IncompressibleMomentumTransportModel
:
    public MomentumTransportModel
    <
        geometricOneField,
        geometricOneField,
        incompressibleMomentumTransportModel,
        TransportModel
    >
{

public:

    typedef geometricOneField alphaField;
    typedef geometricOneField rhoField;
    typedef TransportModel transportModel;


    // Constructors

        //- Construct
        IncompressibleMomentumTransportModel
        (
            const word& type,
            const geometricOneField& alpha,
            const geometricOneField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const TransportModel& transport
        );


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<IncompressibleMomentumTransportModel> New
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            const TransportModel& transportModel
        );


    //- Destructor
    virtual ~IncompressibleMomentumTransportModel()
    {}


    // Member Functions

        //- Return the laminar viscosity
        virtual tmp<volScalarField> nu() const
        {
            return this->transport_.nu();
        }

        //- Return the laminar viscosity on patchi
        virtual tmp<scalarField> nu(const label patchi) const
        {
            return this->transport_.nu(patchi);
        }

        //- Return the effective stress tensor
        virtual tmp<volSymmTensorField> devSigma() const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevSigma(volVectorField& U) const;

        //- Return the effective stress tensor
        virtual tmp<volSymmTensorField> devTau() const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevTau(volVectorField& U) const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevTau
        (
            const volScalarField& rho,
            volVectorField& U
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "IncompressibleMomentumTransportModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
