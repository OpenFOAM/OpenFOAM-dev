/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2023 OpenFOAM Foundation
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
    Foam::segregatedDisplacedSidedPhaseInterface

Description
    Class to represent a certain side of an interface between phases, where one
    phase is considered segregated within the other, and the pair has been
    displaced to some extent by a third phase.

SourceFiles
    segregatedDisplacedSidedPhaseInterface.C

\*---------------------------------------------------------------------------*/

#ifndef segregatedDisplacedSidedPhaseInterface_H
#define segregatedDisplacedSidedPhaseInterface_H

#include "segregatedPhaseInterface.H"
#include "displacedPhaseInterface.H"
#include "sidedPhaseInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class segregatedDisplacedSidedPhaseInterface Declaration
\*---------------------------------------------------------------------------*/

class segregatedDisplacedSidedPhaseInterface
:
    public segregatedPhaseInterface,
    public displacedPhaseInterface,
    public sidedPhaseInterface
{
public:

    //- Runtime type information
    TypeName("segregatedDisplacedSidedPhaseInterface");


    // Constructors

        //- Construct from phases
        segregatedDisplacedSidedPhaseInterface
        (
            const phaseModel& phase,
            const phaseModel& otherPhase,
            const phaseModel& displacing
        );

        //- Construct from fluid and name
        segregatedDisplacedSidedPhaseInterface
        (
            const phaseSystem& fluid,
            const word& name
        );


    //- Destructor
    virtual ~segregatedDisplacedSidedPhaseInterface();


    // Member Functions

        //- Name
        virtual word name() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
