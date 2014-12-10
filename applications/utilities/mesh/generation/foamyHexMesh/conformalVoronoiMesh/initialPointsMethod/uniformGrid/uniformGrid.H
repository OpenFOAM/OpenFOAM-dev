/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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
    Foam::uniformGrid

Description
    Generate a uniform grid of points inside the surfaces to be
    conformed to of the conformalVoronoiMesh

SourceFiles
    uniformGrid.C

\*---------------------------------------------------------------------------*/

#ifndef uniformGrid_H
#define uniformGrid_H

#include "initialPointsMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class uniformGrid Declaration
\*---------------------------------------------------------------------------*/

class uniformGrid
:
    public initialPointsMethod
{

private:

    // Private data

    //- The initial cell spacing
    scalar initialCellSize_;

    //- Should the initial positions be randomised
    Switch randomiseInitialGrid_;

    //- Randomise the initial positions by fraction of the initialCellSize_
    scalar randomPerturbationCoeff_;


public:

    //- Runtime type information
    TypeName("uniformGrid");

    // Constructors

        //- Construct from components
        uniformGrid
        (
            const dictionary& initialPointsDict,
            const Time& runTime,
            Random& rndGen,
            const conformationSurfaces& geometryToConformTo,
            const cellShapeControl& cellShapeControls,
            const autoPtr<backgroundMeshDecomposition>& decomposition
        );


    //- Destructor
    virtual ~uniformGrid()
    {}


    // Member Functions

        //- Return the initial points for the conformalVoronoiMesh
        virtual List<Vb::Point> initialPoints() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
