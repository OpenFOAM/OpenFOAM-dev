/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::processorCyclicGAMGInterfaceField

Description
    GAMG agglomerated processor interface field.

SourceFiles
    processorCyclicGAMGInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef processorCyclicGAMGInterfaceField_H
#define processorCyclicGAMGInterfaceField_H

#include "processorGAMGInterfaceField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class processorCyclicGAMGInterfaceField Declaration
\*---------------------------------------------------------------------------*/

class processorCyclicGAMGInterfaceField
:
    public processorGAMGInterfaceField
{
    // Private Member Functions

        //- Disallow default bitwise copy construction
        processorCyclicGAMGInterfaceField
        (
            const processorCyclicGAMGInterfaceField&
        );

        //- Disallow default bitwise assignment
        void operator=(const processorCyclicGAMGInterfaceField&) = delete;


public:

    //- Runtime type information
    TypeName("processorCyclic");


    // Constructors

        //- Construct from GAMG interface and fine level interface field
        processorCyclicGAMGInterfaceField
        (
            const GAMGInterface& GAMGCp,
            const lduInterfaceField& fineInterface
        );

        //- Construct from GAMG interface and fine level interface field
        processorCyclicGAMGInterfaceField
        (
            const GAMGInterface& GAMGCp,
            const int rank
        );


    // Destructor

        virtual ~processorCyclicGAMGInterfaceField();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
