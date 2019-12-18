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
    Foam::nearWallDistNoSearch

Description
    Distance calculation for cells with face on a wall. Does not search
    anything, just takes normal component of distance.

SourceFiles
    nearWallDistNoSearch.C

\*---------------------------------------------------------------------------*/

#ifndef nearWallDistNoSearch_H
#define nearWallDistNoSearch_H

#include "volFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                           Class nearWallDistNoSearch Declaration
\*---------------------------------------------------------------------------*/

class nearWallDistNoSearch
:
    public volScalarField::Boundary
{
    // Private Data

        //- Reference to mesh
        const fvMesh& mesh_;


    // Private Member Functions

        //- Do all calculations.
        void doAll();


public:

    // Constructors

        //- Construct from components
        nearWallDistNoSearch(const fvMesh& mesh);

        //- Disallow default bitwise copy construction
        nearWallDistNoSearch(const nearWallDistNoSearch&) = delete;


    //- Destructor
    virtual ~nearWallDistNoSearch();


    // Member Functions

        //- Correct for mesh geom/topo changes
        virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const nearWallDistNoSearch&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
