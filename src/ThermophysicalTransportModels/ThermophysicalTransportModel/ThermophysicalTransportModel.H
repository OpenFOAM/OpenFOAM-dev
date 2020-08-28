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
    Foam::ThermophysicalTransportModel

Description
    Templated abstract base class for thermophysical transport models

SourceFiles
    ThermophysicalTransportModel.C

\*---------------------------------------------------------------------------*/

#ifndef ThermophysicalTransportModel_H
#define ThermophysicalTransportModel_H

#include "thermophysicalTransportModel.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class ThermophysicalTransportModel Declaration
\*---------------------------------------------------------------------------*/

template<class MomentumTransportModel, class ThermoModel>
class ThermophysicalTransportModel
:
    public thermophysicalTransportModel
{

public:

    typedef MomentumTransportModel momentumTransportModel;
    typedef ThermoModel thermoModel;
    typedef typename momentumTransportModel::alphaField alphaField;


protected:

    // Protected data

        const momentumTransportModel& momentumTransport_;
        const thermoModel& thermo_;


public:

    // Declare run-time constructor selection table

        declareRunTimeNewSelectionTable
        (
            autoPtr,
            ThermophysicalTransportModel,
            dictionary,
            (
                const momentumTransportModel& momentumTransport,
                const thermoModel& thermo
            ),
            (momentumTransport, thermo)
        );


    // Constructors

        //- Construct
        ThermophysicalTransportModel
        (
            const momentumTransportModel& momentumTransport,
            const thermoModel& thermo
        );

        //- Disallow default bitwise copy construction
        ThermophysicalTransportModel
        (
            const ThermophysicalTransportModel&
        ) = delete;


    // Selectors

        //- Return a reference to the selected thermophysical transport model
        static autoPtr<ThermophysicalTransportModel> New
        (
            const momentumTransportModel& transport,
            const thermoModel& thermo
        );


    //- Destructor
    virtual ~ThermophysicalTransportModel()
    {}


    // Member Functions

        //- Access function to incompressible transport model
        const momentumTransportModel& momentumTransport() const
        {
            return momentumTransport_;
        }

        //- Access function to incompressible transport model
        virtual const thermoModel& thermo() const
        {
            return thermo_;
        }

        //- Access function to phase fraction
        const alphaField& alpha() const
        {
            return momentumTransport_.alpha();
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const ThermophysicalTransportModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ThermophysicalTransportModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
