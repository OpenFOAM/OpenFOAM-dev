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
    Foam::segregatedSidedPhaseInterface

Description
    Class to represent a certain side of an interface between phases, where one
    phase is considered segregated within the other.

SourceFiles
    segregatedSidedPhaseInterface.C

\*---------------------------------------------------------------------------*/

#ifndef segregatedSidedPhaseInterface_H
#define segregatedSidedPhaseInterface_H

#include "segregatedPhaseInterface.H"
#include "sidedPhaseInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class segregatedSidedPhaseInterface Declaration
\*---------------------------------------------------------------------------*/

class segregatedSidedPhaseInterface
:
    public segregatedPhaseInterface,
    public sidedPhaseInterface
{
public:

    //- Runtime type information
    TypeName("segregatedSidedPhaseInterface");


    // Constructors

        //- Construct from phases
        segregatedSidedPhaseInterface
        (
            const phaseModel& phase,
            const phaseModel& otherPhase
        );

        //- Construct from fluid and name
        segregatedSidedPhaseInterface
        (
            const phaseSystem& fluid,
            const word& name
        );


    //- Destructor
    virtual ~segregatedSidedPhaseInterface();


    // Member Functions

        //- Name
        virtual word name() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
