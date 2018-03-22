/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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
    Foam::AnisothermalPhaseModel

Description
    Class which represents a phase for which the temperature (strictly energy)
    varies. Returns the energy equation and corrects the thermodynamic model.

SourceFiles
    AnisothermalPhaseModel.C

\*---------------------------------------------------------------------------*/

#ifndef AnisothermalPhaseModel_H
#define AnisothermalPhaseModel_H

#include "phaseModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class AnisothermalPhaseModel Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseModel>
class AnisothermalPhaseModel
:
    public BasePhaseModel
{
    // Private member functions

        //- Optionally filter the pressure work term as the phase-fraction -> 0
        tmp<volScalarField> filterPressureWork
        (
            const tmp<volScalarField>& pressureWork
        ) const;


public:

    // Constructors

        AnisothermalPhaseModel
        (
            const phaseSystem& fluid,
            const word& phaseName,
            const label index
        );


    //- Destructor
    virtual ~AnisothermalPhaseModel();


    // Member Functions

        //- Correct the thermodynamics
        virtual void correctThermo();

        //- Return whether the phase is isothermal
        virtual bool isothermal() const;

        //- Return the enthalpy equation
        virtual tmp<fvScalarMatrix> heEqn();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "AnisothermalPhaseModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
