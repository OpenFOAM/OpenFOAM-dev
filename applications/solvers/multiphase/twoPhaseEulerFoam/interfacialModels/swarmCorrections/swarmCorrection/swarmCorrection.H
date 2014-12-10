/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 OpenFOAM Foundation
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

SourceFiles
    swarmCorrection.C
    newSwarmCorrection.C

\*---------------------------------------------------------------------------*/

#ifndef swarmCorrection_H
#define swarmCorrection_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class phasePair;

/*---------------------------------------------------------------------------*\
                           Class swarmCorrection Declaration
\*---------------------------------------------------------------------------*/

class swarmCorrection
{
protected:

    // Protected data

        //- Phase pair
        const phasePair& pair_;


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
            const phasePair& pair
        ),
        (dict, pair)
    );


    // Constructors

        //- Construct from a dictionary and a phase pair
        swarmCorrection
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~swarmCorrection();


    // Selectors

        static autoPtr<swarmCorrection> New
        (
            const dictionary& dict,
            const phasePair& pair
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
