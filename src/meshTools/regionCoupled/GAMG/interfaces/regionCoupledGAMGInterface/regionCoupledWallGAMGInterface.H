/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::regionCoupledWallGAMGInterface

Description
    GAMG agglomerated coupled region interface.

SourceFiles
    regionCoupledWallGAMGInterface.C

\*---------------------------------------------------------------------------*/

#ifndef regionCoupledWallGAMGInterface_H
#define regionCoupledWallGAMGInterface_H

#include "regionCoupledBaseGAMGInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class regionCoupledWallGAMGInterface Declaration
\*---------------------------------------------------------------------------*/

class regionCoupledWallGAMGInterface
:
    public regionCoupledBaseGAMGInterface
{
public:

    //- Runtime type information
    TypeName("regionCoupledWall");


    // Constructors

        //- Construct from fine level interface,
        //  local and neighbour restrict addressing
        regionCoupledWallGAMGInterface
        (
            const label index,
            const lduInterfacePtrsList& coarseInterfaces,
            const lduInterface& fineInterface,
            const labelField& restrictAddressing,
            const labelField& neighbourRestrictAddressing,
            const label fineLevelIndex,
            const label coarseComm
        );

        //- Disallow default bitwise copy construction
        regionCoupledWallGAMGInterface
        (
            const regionCoupledWallGAMGInterface&
        ) = delete;


    //- Destructor
    virtual ~regionCoupledWallGAMGInterface();


    // Member Functions

        // I/O

            //- Write to stream
            virtual void write(Ostream&) const
            {
                // TBD. How to serialise the AMI such that we can stream
                // regionCoupledWallGAMGInterface.
                NotImplemented;
            }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const regionCoupledWallGAMGInterface&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
