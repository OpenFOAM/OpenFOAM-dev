/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2023 OpenFOAM Foundation
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
    Foam::swarmCorrection

Description
    Drag model modifier to take into account large fractions (swarms) of the
    dispersed phase. Often not appropriate as these effects are a fundamental
    part of the drag model formulation (e.g. Ergun). Use with care.

SourceFiles
    swarmCorrection.C
    swarmCorrectionNew.C

\*---------------------------------------------------------------------------*/

#ifndef swarmCorrection_H
#define swarmCorrection_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "dispersedPhaseInterface.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class swarmCorrection Declaration
\*---------------------------------------------------------------------------*/

class swarmCorrection
{
protected:

    // Protected data

        //- Interface
        const dispersedPhaseInterface interface_;


public:

    //- Runtime type information
    TypeName("swarmCorrection");


    // Declare runtime construction
    declareRunTimeSelectionTable
    (
        autoPtr,
        swarmCorrection,
        dictionary,
        (
            const dictionary& dict,
            const phaseInterface& interface
        ),
        (dict, interface)
    );


    // Constructors

        //- Construct from a dictionary and an interface
        swarmCorrection
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    //- Destructor
    virtual ~swarmCorrection();


    // Selectors

        static autoPtr<swarmCorrection> New
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    // Member Functions

        //- Swarm correction coefficient
        virtual tmp<volScalarField> Cs() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
