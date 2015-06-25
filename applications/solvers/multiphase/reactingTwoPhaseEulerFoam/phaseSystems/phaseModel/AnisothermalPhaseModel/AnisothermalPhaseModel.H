/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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
                           Class phaseModel Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseModel>
class AnisothermalPhaseModel
:
    public BasePhaseModel
{
    // Private data

        //- Dilatation
        volScalarField D_;

        //- Kinetic energy
        volScalarField K_;


public:

    // Constructors

        AnisothermalPhaseModel(const phaseSystem& fluid, const word& phaseName);


    //- Destructor
    virtual ~AnisothermalPhaseModel();


    // Member Functions

        //- Correct the kinematics
        virtual void correctKinematics();

        //- Correct the thermodynamics
        virtual void correctThermo();

        //- Return the enthalpy equation
        virtual tmp<fvScalarMatrix> heEqn();


        // Compressibility (variable density)

            //- Return true if the phase is compressible otherwise false
            virtual bool compressible() const;

            //- Phase dilatation rate ((alpha/rho)*Drho/Dt)
            virtual tmp<volScalarField> D() const;

            //- Set phase dilatation rate ((alpha/rho)*Drho/Dt)
            virtual void D(const volScalarField& D);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "AnisothermalPhaseModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
