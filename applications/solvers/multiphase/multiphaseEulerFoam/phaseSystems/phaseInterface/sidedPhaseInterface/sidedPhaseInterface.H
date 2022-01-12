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
    Foam::sidedPhaseInterface

Description
    Class to represent a certain side of an interface between phases.

SourceFiles
    sidedPhaseInterface.C

\*---------------------------------------------------------------------------*/

#ifndef sidedPhaseInterface_H
#define sidedPhaseInterface_H

#include "phaseInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class sidedPhaseInterface Declaration
\*---------------------------------------------------------------------------*/

class sidedPhaseInterface
:
    virtual public phaseInterface
{
private:

    // Private Data

        //- Phase
        const phaseModel& phase_;


public:

    //- Runtime type information
    TypeName("sidedPhaseInterface");


    // Constructors

        //- Construct from phases
        sidedPhaseInterface
        (
            const phaseModel& phase,
            const phaseModel& otherPhase
        );

        //- Construct from fluid and name
        sidedPhaseInterface
        (
            const phaseSystem& fluid,
            const word& name
        );


    //- Destructor
    virtual ~sidedPhaseInterface();


    // Static Member Functions

        //- Return the separator that delimits this interface's name
        static word separator()
        {
            return "inThe";
        }


    // Member Functions

        //- Name
        virtual word name() const;


        // Access

            //- Phase on the this side of the interface
            const phaseModel& phase() const;

            //- Phase on the other side of the interface
            const phaseModel& otherPhase() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
