/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    Foam::IsothermalSolidPhaseModel

Description
    Class which represents a solid phase for which the temperature (strictly
    energy) remains constant. Returns an empty energy equation and does nothing
    when correctThermo is called.

SourceFiles
    IsothermalSolidPhaseModel.C

\*---------------------------------------------------------------------------*/

#ifndef IsothermalSolidPhaseModel_H
#define IsothermalSolidPhaseModel_H

#include "phaseModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class IsothermalSolidPhaseModel Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseModel>
class IsothermalSolidPhaseModel
:
    public BasePhaseModel
{
public:

    // Constructors

        IsothermalSolidPhaseModel
        (
            const phaseSystem& fluid,
            const word& phaseName,
            const bool referencePhase,
            const label index
        );


    //- Destructor
    virtual ~IsothermalSolidPhaseModel();


    // Member Functions

        //- Correct the thermodynamics
        virtual void correctThermo();

        //- Return whether the phase is isothermal
        virtual bool isothermal() const;

        //- Return the effective thermal conductivity on a patch
        virtual tmp<scalarField> kappaEff(const label patchi) const;

        //- Return the enthalpy equation
        virtual tmp<fvScalarMatrix> heEqn();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "IsothermalSolidPhaseModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
