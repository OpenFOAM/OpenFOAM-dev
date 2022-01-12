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
    Foam::segregatedPhaseInterface

Description
    Class to represent a interface between phases where the two phases are
    considered to be segregated; that is, separated by a geometrically complex
    interface for which dispersed representations are inappropriate.

SourceFiles
    segregatedPhaseInterface.C

\*---------------------------------------------------------------------------*/

#ifndef segregatedPhaseInterface_H
#define segregatedPhaseInterface_H

#include "phaseInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class segregatedPhaseInterface Declaration
\*---------------------------------------------------------------------------*/

class segregatedPhaseInterface
:
    virtual public phaseInterface
{
public:

    //- Runtime type information
    TypeName("segregatedPhaseInterface");


    // Constructors

        //- Construct from phases
        segregatedPhaseInterface
        (
            const phaseModel& phase1,
            const phaseModel& phase2
        );

        //- Construct from fluid and name
        segregatedPhaseInterface
        (
            const phaseSystem& fluid,
            const word& name
        );


    //- Destructor
    virtual ~segregatedPhaseInterface();


    // Static Member Functions

        //- Return the separator that delimits this interface's name
        static word separator()
        {
            return "segregatedWith";
        }


    // Member Functions

        //- Name
        virtual word name() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
