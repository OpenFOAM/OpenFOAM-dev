/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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
    Foam::diameterModels::driftModels::isothermal

Description
    Drift induced by interfacial phaseChange.

SourceFiles
    isothermal.C

\*---------------------------------------------------------------------------*/

#ifndef phaseChange_H
#define phaseChange_H

#include "driftModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace driftModels
{

/*---------------------------------------------------------------------------*\
                     Class phaseChange Declaration
\*---------------------------------------------------------------------------*/

class phaseChange
:
    public driftModel
{
    // Private data

        //- Names of unordered phasePairs between which phaseChange occurs, i.e.
        //  "(gasIAndLiquid gasIIAndLiquid)"
        List<word> pairNames_;

        //- Total mass transfer rate due to phaseChange
        volScalarField iDmdt_;

        //- Total number concentration
        volScalarField N_;

public:

    //- Runtime type information
    TypeName("phaseChange");

    // Constructor

        phaseChange
        (
            const populationBalanceModel& popBal,
            const dictionary& dict
        );


    //- Destructor
    virtual ~phaseChange()
    {}


    // Member Functions

        //- Correct diameter independent expressions
        virtual void correct();

        //- Add to driftRate
        virtual void addToDriftRate
        (
            volScalarField& driftRate,
            const label i
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace driftModels
} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
