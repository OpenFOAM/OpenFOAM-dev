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
    Foam::cyclicAMILduInterface

Description
    An abstract base class for cyclic AMI coupled interfaces

SourceFiles
    cyclicAMILduInterface.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicAMILduInterface_H
#define cyclicAMILduInterface_H

#include "primitiveFieldsFwd.H"
#include "AMIInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class cyclicAMILduInterface Declaration
\*---------------------------------------------------------------------------*/

class cyclicAMILduInterface
{

public:

    //- Runtime type information
    TypeName("cyclicAMILduInterface");


    // Constructors

        //- Construct null
        cyclicAMILduInterface()
        {}


    //- Destructor
    virtual ~cyclicAMILduInterface();


    // Member Functions

        // Access

            //- Return neighbour
            virtual label nbrPatchID() const = 0;

            //- Does this side own the interface?
            virtual bool owner() const = 0;

            //- Return processor number
            virtual const cyclicAMILduInterface& nbrPatch() const = 0;

            //- Return a reference to the AMI interpolators
            virtual const PtrList<AMIInterpolation>& AMIs() const = 0;

            // Return a reference to the AMI transformations
            virtual const List<transformer>&
                AMITransforms() const = 0;

            //- Return transformation between the coupled patches
            virtual const transformer& transform() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
