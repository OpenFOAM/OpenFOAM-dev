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
    Foam::dispersedDisplacedPhaseInterface

Description
    Class to represent a interface between phases where one phase is considered
    dispersed within the other, and the interface has been displaced to some
    extent by a third phase.

SourceFiles
    dispersedDisplacedPhaseInterface.C

\*---------------------------------------------------------------------------*/

#ifndef dispersedDisplacedPhaseInterface_H
#define dispersedDisplacedPhaseInterface_H

#include "dispersedPhaseInterface.H"
#include "displacedPhaseInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class dispersedDisplacedPhaseInterface Declaration
\*---------------------------------------------------------------------------*/

class dispersedDisplacedPhaseInterface
:
    public dispersedPhaseInterface,
    public displacedPhaseInterface
{
public:

    //- Runtime type information
    TypeName("dispersedDisplacedPhaseInterface");


    // Constructors

        //- Construct from phases
        dispersedDisplacedPhaseInterface
        (
            const phaseModel& dispersed,
            const phaseModel& continuous,
            const phaseModel& displacing
        );

        //- Construct from fluid and name
        dispersedDisplacedPhaseInterface
        (
            const phaseSystem& fluid,
            const word& name
        );


    //- Destructor
    virtual ~dispersedDisplacedPhaseInterface();


    // Member Functions

        //- Name
        virtual word name() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
