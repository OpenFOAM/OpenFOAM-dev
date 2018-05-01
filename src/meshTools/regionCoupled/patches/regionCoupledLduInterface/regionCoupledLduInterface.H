/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::regionCoupledLduInterface

Description
    An abstract base class for region coupled interfaces

SourceFiles
    regionCoupledLduInterface.C

\*---------------------------------------------------------------------------*/

#ifndef regionCoupledLduInterface_H
#define regionCoupledLduInterface_H

#include "primitiveFieldsFwd.H"
#include "AMIPatchToPatchInterpolation.H"
#include "polyMesh.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class regionCoupledLduInterface Declaration
\*---------------------------------------------------------------------------*/

class regionCoupledLduInterface
{

public:

    //- Runtime type information
    TypeName("regionCoupledLduInterface");


    // Constructors

        //- Construct null
        regionCoupledLduInterface()
        {}


    //- Destructor
    virtual ~regionCoupledLduInterface();


    // Member Functions

        // Access

            //- Return neighbour
            virtual label neighbPatchID() const = 0;

            //- Is it owner?
            virtual bool owner() const = 0;

            //- Return neighb regionCoupledLduInterface
            virtual const regionCoupledLduInterface& neighbPatch() const = 0;

            //- Return AMI
            // virtual const AMIPatchToPatchInterpolation& AMI() const = 0;

            //- Return nbrMesh
            virtual const polyMesh& nbrMesh() const = 0;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
