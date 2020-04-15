/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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
    Foam::LESThermophysicalTransportModel

Description
    Templated abstract base class for LES thermophysical transport models

SourceFiles
    LESThermophysicalTransportModel.C

\*---------------------------------------------------------------------------*/

#ifndef LESThermophysicalTransportModel_H
#define LESThermophysicalTransportModel_H

#include "ThermophysicalTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class LESThermophysicalTransportModel Declaration
\*---------------------------------------------------------------------------*/

template<class BasicThermophysicalTransportModel>
class LESThermophysicalTransportModel
:
    public BasicThermophysicalTransportModel
{

protected:

    // Protected data

        //- LES coefficients dictionary
        dictionary LESDict_;

        //- Flag to print the model coeffs at run-time
        Switch printCoeffs_;

        //- Model coefficients dictionary
        dictionary coeffDict_;


    // Protected Member Functions

        //- Print model coefficients
        virtual void printCoeffs(const word& type);


public:

    typedef typename BasicThermophysicalTransportModel::alphaField
        alphaField;

    typedef typename BasicThermophysicalTransportModel::momentumTransportModel
        momentumTransportModel;


    //- Runtime type information
    TypeName("LES");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            LESThermophysicalTransportModel,
            dictionary,
            (
                const momentumTransportModel& momentumTransport
            ),
            (momentumTransport)
        );


    // Constructors

        //- Construct from components
        LESThermophysicalTransportModel
        (
            const word& type,
            const momentumTransportModel& momentumTransport
        );

        //- Disallow default bitwise copy construction
        LESThermophysicalTransportModel
        (
            const LESThermophysicalTransportModel&
        ) = delete;


    // Selectors

        //- Return a reference to the selected LES model
        static autoPtr<LESThermophysicalTransportModel> New
        (
            const momentumTransportModel& momentumTransport
        );


    //- Destructor
    virtual ~LESThermophysicalTransportModel()
    {}


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read();

        //- Const access to the coefficients dictionary
        virtual const dictionary& coeffDict() const
        {
            return coeffDict_;
        }

        //- Thermal diffusivity for temperature of mixture [W/m/K]
        virtual tmp<volScalarField> kappa() const
        {
            return this->thermo().kappa();
        }

        //- Thermal diffusivity for temperature of mixture
        //  for patch [W/m/K]
        virtual tmp<scalarField> kappa(const label patchi) const
        {
            return this->thermo().kappa(patchi);
        }

        //- Effective thermal turbulent diffusivity for temperature
        //  of mixture [W/m/K]
        virtual tmp<volScalarField> kappaEff() const
        {
            return volScalarField::New
            (
                "kappaEff",
                this->thermo().kappa()
            );
        }

        //- Effective thermal turbulent diffusivity for temperature
        //  of mixture for patch [W/m/K]
        virtual tmp<scalarField> kappaEff(const label patchi) const
        {
            return kappa(patchi);
        }

        //- Effective thermal turbulent diffusivity of mixture [kg/m/s]
        virtual tmp<volScalarField> alphaEff() const
        {
            return volScalarField::New
            (
                "alphaEff",
                this->thermo().alphahe()
            );
        }

        //- Effective thermal turbulent diffusivity of mixture
        //  for patch [kg/m/s]
        virtual tmp<scalarField> alphaEff(const label patchi) const
        {
            return this->thermo().alphahe(patchi);
        }

        //- Correct the LES transport
        virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const LESThermophysicalTransportModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "LESThermophysicalTransportModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
