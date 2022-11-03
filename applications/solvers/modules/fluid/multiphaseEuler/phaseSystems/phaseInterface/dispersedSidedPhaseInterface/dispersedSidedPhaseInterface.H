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
    Foam::dispersedSidedPhaseInterface

Description
    Class to represent a certain side of an interface between phases where one
    phase is considered dispersed within the other.

SourceFiles
    dispersedSidedPhaseInterface.C

\*---------------------------------------------------------------------------*/

#ifndef dispersedSidedPhaseInterface_H
#define dispersedSidedPhaseInterface_H

#include "dispersedPhaseInterface.H"
#include "sidedPhaseInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class dispersedSidedPhaseInterface Declaration
\*---------------------------------------------------------------------------*/

class dispersedSidedPhaseInterface
:
    public dispersedPhaseInterface,
    public sidedPhaseInterface
{
public:

    //- Runtime type information
    TypeName("dispersedSidedPhaseInterface");


    // Constructors

        //- Construct from phases
        dispersedSidedPhaseInterface
        (
            const phaseModel& dispersed,
            const phaseModel& continuous,
            const phaseModel& phase
        );

        //- Construct from fluid and name
        dispersedSidedPhaseInterface
        (
            const phaseSystem& fluid,
            const word& name
        );


    //- Destructor
    virtual ~dispersedSidedPhaseInterface();


    // Member Functions

        //- Name
        virtual word name() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
