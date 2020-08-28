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
    Foam::RASThermophysicalTransportModel

Description
    Templated abstract base class for RAS thermophysical transport models

SourceFiles
    RASThermophysicalTransportModel.C

\*---------------------------------------------------------------------------*/

#ifndef RASThermophysicalTransportModel_H
#define RASThermophysicalTransportModel_H

#include "ThermophysicalTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class RASThermophysicalTransportModel Declaration
\*---------------------------------------------------------------------------*/

template<class BasicThermophysicalTransportModel>
class RASThermophysicalTransportModel
:
    public BasicThermophysicalTransportModel
{

protected:

    // Protected data

        //- RAS coefficients dictionary
        dictionary RASDict_;

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

    typedef typename BasicThermophysicalTransportModel::thermoModel
        thermoModel;


    //- Runtime type information
    TypeName("RAS");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            RASThermophysicalTransportModel,
            dictionary,
            (
                const momentumTransportModel& momentumTransport,
                const thermoModel& thermo
            ),
            (momentumTransport, thermo)
        );


    // Constructors

        //- Construct from components
        RASThermophysicalTransportModel
        (
            const word& type,
            const momentumTransportModel& momentumTransport,
            const thermoModel& thermo
        );

        //- Disallow default bitwise copy construction
        RASThermophysicalTransportModel
        (
            const RASThermophysicalTransportModel&
        ) = delete;


    // Selectors

        //- Return a reference to the selected RAS model
        static autoPtr<RASThermophysicalTransportModel> New
        (
            const momentumTransportModel& momentumTransport,
            const thermoModel& thermo
        );


    //- Destructor
    virtual ~RASThermophysicalTransportModel()
    {}


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read();

        //- Const access to the coefficients dictionary
        virtual const dictionary& coeffDict() const
        {
            return coeffDict_;
        }

        //- Effective thermal turbulent diffusivity for temperature
        //  of mixture [W/m/K]
        virtual tmp<volScalarField> kappaEff() const = 0;

        //- Effective thermal turbulent diffusivity for temperature
        //  of mixture for patch [W/m/K]
        virtual tmp<scalarField> kappaEff(const label patchi) const = 0;

        //- Effective thermal turbulent diffusivity of mixture [kg/m/s]
        virtual tmp<volScalarField> alphaEff() const = 0;

        //- Effective thermal turbulent diffusivity of mixture
        //  for patch [kg/m/s]
        virtual tmp<scalarField> alphaEff(const label patchi) const = 0;

        //- Correct the RAS transport
        virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const RASThermophysicalTransportModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "RASThermophysicalTransportModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
