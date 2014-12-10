/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::reflectionVectors

Description
    Container for reflection vectors (= unit normal of nearest wall)
SourceFiles
    reflectionVectors.C

\*---------------------------------------------------------------------------*/

#ifndef reflectionVectors_H
#define reflectionVectors_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class reflectionVectors Declaration
\*---------------------------------------------------------------------------*/

class reflectionVectors
{

protected:

    // Protected data

        volVectorField n_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        reflectionVectors(const reflectionVectors&);

        //- Disallow default bitwise assignment
        void operator=(const reflectionVectors&);


public:

    // Constructors

        //- Construct from mesh
        reflectionVectors(const fvMesh& mesh);


    // Member Functions

        const volVectorField& n() const
        {
            return n_;
        }


        //- Correct for mesh geom/topo changes
        void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
