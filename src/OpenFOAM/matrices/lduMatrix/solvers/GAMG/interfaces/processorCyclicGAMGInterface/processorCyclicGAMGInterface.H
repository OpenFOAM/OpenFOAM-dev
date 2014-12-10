/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::processorCyclicGAMGInterface

Description
    GAMG agglomerated processor interface.

SourceFiles
    processorCyclicGAMGInterface.C
    processorCyclicGAMGInterfaceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef processorCyclicGAMGInterface_H
#define processorCyclicGAMGInterface_H

#include "processorGAMGInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class processorCyclicGAMGInterface Declaration
\*---------------------------------------------------------------------------*/

class processorCyclicGAMGInterface
:
    public processorGAMGInterface
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        processorCyclicGAMGInterface(const processorCyclicGAMGInterface&);

        //- Disallow default bitwise assignment
        void operator=(const processorCyclicGAMGInterface&);


public:

    //- Runtime type information
    TypeName("processorCyclic");

    // Constructors

        //- Construct from fine-level interface,
        //  local and neighbour restrict addressing
        processorCyclicGAMGInterface
        (
            const label index,
            const lduInterfacePtrsList& coarseInterfaces,
            const lduInterface& fineInterface,
            const labelField& restrictAddressing,
            const labelField& neighbourRestrictAddressing,
            const label fineLevelIndex,
            const label coarseComm
        );

        //- Construct from Istream
        processorCyclicGAMGInterface
        (
            const label index,
            const lduInterfacePtrsList& coarseInterfaces,
            Istream& is
        );


    // Destructor

        virtual ~processorCyclicGAMGInterface();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
