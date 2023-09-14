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
    Foam::dispersedDisplacedSidedPhaseInterface

Description
    Class to represent a certain side of an interface between phases, where one
    phase is considered dispersed within the other, and the interface has been
    displaced to some extent by a third phase.

SourceFiles
    dispersedDisplacedSidedPhaseInterface.C

\*---------------------------------------------------------------------------*/

#ifndef dispersedDisplacedSidedPhaseInterface_H
#define dispersedDisplacedSidedPhaseInterface_H

#include "dispersedPhaseInterface.H"
#include "displacedPhaseInterface.H"
#include "sidedPhaseInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class dispersedDisplacedSidedPhaseInterface Declaration
\*---------------------------------------------------------------------------*/

class dispersedDisplacedSidedPhaseInterface
:
    public dispersedPhaseInterface,
    public displacedPhaseInterface,
    public sidedPhaseInterface
{
public:

    //- Runtime type information
    TypeName("dispersedDisplacedSidedPhaseInterface");


    // Constructors

        //- Construct from phases
        dispersedDisplacedSidedPhaseInterface
        (
            const phaseModel& dispersed,
            const phaseModel& continuous,
            const phaseModel& displacing,
            const phaseModel& phase
        );

        //- Construct from fluid and name
        dispersedDisplacedSidedPhaseInterface
        (
            const phaseSystem& fluid,
            const word& name
        );


    //- Destructor
    virtual ~dispersedDisplacedSidedPhaseInterface();


    // Member Functions

        //- Name
        virtual word name() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
