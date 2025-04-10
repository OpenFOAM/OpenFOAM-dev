/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2024 OpenFOAM Foundation
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
    Foam::laminarThermophysicalTransportModel

Description
    Templated abstract base class for laminar thermophysical transport models

SourceFiles
    laminarThermophysicalTransportModel.C

\*---------------------------------------------------------------------------*/

#ifndef laminarThermophysicalTransportModel_H
#define laminarThermophysicalTransportModel_H

#include "ThermophysicalTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class laminarThermophysicalTransportModel Declaration
\*---------------------------------------------------------------------------*/

template<class BasicThermophysicalTransportModel>
class laminarThermophysicalTransportModel
:
    public BasicThermophysicalTransportModel
{

protected:

    // Protected member functions

        //- Const access to the coefficients dictionary
        const dictionary& coeffDict() const;

        //- Effective thermal turbulent diffusivity of mixture [kg/m/s]
        //  Used for the implicit energy correction on the temperature laplacian
        tmp<volScalarField> alphaEff() const
        {
            return volScalarField::New
            (
                "alphaEff",
                this->thermo().kappa()/this->thermo().Cpv()
            );
        }


public:

    typedef typename BasicThermophysicalTransportModel::alphaField
        alphaField;

    typedef typename BasicThermophysicalTransportModel::momentumTransportModel
        momentumTransportModel;

    typedef typename BasicThermophysicalTransportModel::thermoModel
        thermoModel;


    //- Runtime type information
    TypeName("laminar");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            laminarThermophysicalTransportModel,
            dictionary,
            (
                const momentumTransportModel& momentumTransport,
                const thermoModel& thermo
            ),
            (momentumTransport, thermo)
        );


    // Constructors

        //- Construct from components
        laminarThermophysicalTransportModel
        (
            const word& type,
            const momentumTransportModel& momentumTransport,
            const thermoModel& thermo
        );

        //- Disallow default bitwise copy construction
        laminarThermophysicalTransportModel
        (
            const laminarThermophysicalTransportModel&
        ) = delete;


    // Selectors

        //- Return a reference to the selected laminar model
        static autoPtr<laminarThermophysicalTransportModel> New
        (
            const momentumTransportModel& momentumTransport,
            const thermoModel& thermo
        );


    //- Destructor
    virtual ~laminarThermophysicalTransportModel()
    {}


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read();

        //- Effective thermal turbulent conductivity
        //  of mixture [W/m/K]
        virtual tmp<volScalarField> kappaEff() const
        {
            return volScalarField::New
            (
                "kappaEff",
                this->thermo().kappa()
            );
        }

        //- Effective thermal turbulent conductivity
        //  of mixture for patch [W/m/K]
        virtual tmp<scalarField> kappaEff(const label patchi) const
        {
            return this->thermo().kappa().boundaryField()[patchi];
        }

        //- Effective thermal turbulent diffusivity
        //  of mixture for a patch [kg/m/s]
        virtual tmp<scalarField> alphaEff(const label patchi) const
        {
            return
                this->thermo().kappa().boundaryField()[patchi]
               /this->thermo().Cpv().boundaryField()[patchi];
        }

        //- Effective mass diffusion coefficient
        //  for a given specie mass-fraction [kg/m/s]
        virtual tmp<volScalarField> DEff(const volScalarField& Yi) const = 0;

        //- Effective mass diffusion coefficient
        //  for a given specie mass-fraction for patch [kg/m/s]
        virtual tmp<scalarField> DEff
        (
            const volScalarField& Yi,
            const label patchi
        ) const = 0;

        //- Predict the laminar transport coefficients if possible
        //  without solving thermophysical transport model equations
        virtual void predict();

        //- Solve the thermophysical transport model equations
        //  and correct the laminar transport coefficients
        virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const laminarThermophysicalTransportModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "laminarThermophysicalTransportModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
