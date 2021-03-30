/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2021 OpenFOAM Foundation
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
    Foam::ThermoPhaseModel

Description
    Class which represents a phase with a thermodynamic model. Provides access
    to the thermodynamic variables. Note that the thermo model itself is not
    returned as this class could be substituted in the hierarchy for one which
    mirrors the functionality, but does not include a thermo model; an
    incompressible phase model, for example.

SourceFiles
    ThermoPhaseModel.C

\*---------------------------------------------------------------------------*/

#ifndef ThermoPhaseModel_H
#define ThermoPhaseModel_H

#include "phaseModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class rhoThermo;

/*---------------------------------------------------------------------------*\
                      Class ThermoPhaseModel Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseModel, class ThermoModel>
class ThermoPhaseModel
:
    public BasePhaseModel,
    public dynamicTransportModel
{
protected:

    // Protected data

        //- Thermophysical model
        autoPtr<ThermoModel> thermo_;


public:

    typedef ThermoModel thermoModel;


    // Constructors

        ThermoPhaseModel
        (
            const phaseSystem& fluid,
            const word& phaseName,
            const bool referencePhase,
            const label index
        );


    //- Destructor
    virtual ~ThermoPhaseModel();


    // Member Functions

        // Thermo

            //- Return whether the phase is incompressible
            virtual bool incompressible() const;

            //- Return whether the phase is constant density
            virtual bool isochoric() const;

            //- Return the thermophysical model
            virtual const rhoThermo& thermo() const;

            //- Access the thermophysical model
            virtual rhoThermo& thermoRef();

            //- Return the density field
            virtual tmp<volScalarField> rho() const;


        // Transport

            //- Return the laminar dynamic viscosity
            virtual tmp<volScalarField> mu() const;

            //- Return the laminar dynamic viscosity on a patch
            virtual tmp<scalarField> mu(const label patchi) const;

            //- Return the laminar kinematic viscosity
            virtual tmp<volScalarField> nu() const;

            //- Return the laminar kinematic viscosity on a patch
            virtual tmp<scalarField> nu(const label patchi) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ThermoPhaseModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
