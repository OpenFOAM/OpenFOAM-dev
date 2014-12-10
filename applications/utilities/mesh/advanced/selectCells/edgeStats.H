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
    Foam::edgeStats

Description
    Helper class to calculate minimum edge length on mesh.

SourceFiles
    edgeStats.C

\*---------------------------------------------------------------------------*/

#ifndef edgeStats_H
#define edgeStats_H

#include "direction.H"
#include "scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class Ostream;
class twoDPointCorrector;

/*---------------------------------------------------------------------------*\
                           Class edgeStats Declaration
\*---------------------------------------------------------------------------*/

class edgeStats
{
    // Private data

        //- Reference to mesh.
        const polyMesh& mesh_;

        //- Component (0,1,2) of normal direction or 3 if 3D case.
        direction normalDir_;

    // Private Member Functions

        //- If 2d get component of normal dir.
        direction getNormalDir(const twoDPointCorrector*) const;


        //- Disallow default bitwise copy construct
        edgeStats(const edgeStats&);

        //- Disallow default bitwise assignment
        void operator=(const edgeStats&);


public:

    // Static data members

        // Max (cos of) angle for edges to be considered aligned with axis.
        static const scalar edgeTol_;


    // Constructors

        //- Construct from mesh
        edgeStats(const polyMesh& mesh);

        //- Construct from mesh and corrector
        edgeStats(const polyMesh& mesh, const twoDPointCorrector* );


    // Member Functions

        //- Calculate minimum edge length and print
        scalar minLen(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
