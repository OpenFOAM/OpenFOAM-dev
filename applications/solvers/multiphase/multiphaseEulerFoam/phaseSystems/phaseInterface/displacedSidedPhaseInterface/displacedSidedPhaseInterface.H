/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2022 OpenFOAM Foundation
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
    Foam::displacedSidedPhaseInterface

Description
    Class to represent a certain side of an interface between phases, where the
    interface has been displaced to some extent by a third phase.

SourceFiles
    displacedSidedPhaseInterface.C

\*---------------------------------------------------------------------------*/

#ifndef displacedSidedPhaseInterface_H
#define displacedSidedPhaseInterface_H

#include "segregatedPhaseInterface.H"
#include "displacedPhaseInterface.H"
#include "sidedPhaseInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class displacedSidedPhaseInterface Declaration
\*---------------------------------------------------------------------------*/

class displacedSidedPhaseInterface
:
    public displacedPhaseInterface,
    public sidedPhaseInterface
{
public:

    //- Runtime type information
    TypeName("displacedSidedPhaseInterface");


    // Constructors

        //- Construct from phases
        displacedSidedPhaseInterface
        (
            const phaseModel& phase,
            const phaseModel& otherPhase,
            const phaseModel& displacing
        );

        //- Construct from fluid and name
        displacedSidedPhaseInterface
        (
            const phaseSystem& fluid,
            const word& name
        );


    //- Destructor
    virtual ~displacedSidedPhaseInterface();


    // Member Functions

        //- Name
        virtual word name() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
