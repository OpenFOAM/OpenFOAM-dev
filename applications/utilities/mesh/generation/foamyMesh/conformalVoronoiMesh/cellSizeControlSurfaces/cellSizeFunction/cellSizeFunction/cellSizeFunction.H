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
    Foam::cellSizeFunction

Description
    Abstract base class for specifying target cell sizes

SourceFiles
    cellSizeFunction.C

\*---------------------------------------------------------------------------*/

#ifndef cellSizeFunction_H
#define cellSizeFunction_H

#include "point.H"
#include "conformalVoronoiMesh.H"
#include "searchableSurface.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "surfaceCellSizeFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class cellSizeFunction Declaration
\*---------------------------------------------------------------------------*/

class cellSizeFunction
:
    public dictionary
{

public:

    //- Surface side mode
    enum sideMode
    {
        smInside,         // Control inside the surface
        smOutside,        // Control outside the surface
        rmBothsides       // Control on both sides of a surface
    };

    //- Runtime type information
    TypeName("cellSizeFunction");


protected:

    // Static data

        //- Point closeness tolerance to a surface where the function "snaps" to
        //  including the surface
        static scalar snapToSurfaceTol_;


    // Protected data

        //- Reference to the searchableSurface that cellSizeFunction
        //  relates to
        const searchableSurface& surface_;

        //- Cell size at the surface
        scalarField surfaceCellSize_;

        autoPtr<surfaceCellSizeFunction> surfaceCellSizeFunction_;

        //- Method details dictionary
        dictionary coeffsDict_;

        const scalar& defaultCellSize_;

        //- Index of the region of the surface that this cell size function
        //  applies to
        const labelList regionIndices_;

        //- Mode of size specification, i.e. inside, outside or bothSides
        sideMode sideMode_;

        label priority_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        cellSizeFunction(const cellSizeFunction&);

        //- Disallow default bitwise assignment
        void operator=(const cellSizeFunction&);


public:

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            cellSizeFunction,
            dictionary,
            (
                const dictionary& cellSizeFunctionDict,
                const searchableSurface& surface,
                const scalar& defaultCellSize,
                const labelList regionIndices
            ),
            (cellSizeFunctionDict, surface, defaultCellSize, regionIndices)
        );


    // Constructors

        //- Construct from components
        cellSizeFunction
        (
            const word& type,
            const dictionary& cellSizeFunctionDict,
            const searchableSurface& surface,
            const scalar& defaultCellSize,
            const labelList regionIndices
        );


    // Selectors

        //- Return a reference to the selected cellSizeFunction
        static autoPtr<cellSizeFunction> New
        (
            const dictionary& cellSizeFunctionDict,
            const searchableSurface& surface,
            const scalar& defaultCellSize,
            const labelList regionIndices
        );


    //- Destructor
    virtual ~cellSizeFunction();


    // Member Functions

        //- Const access to the details dictionary
        inline const dictionary& coeffsDict() const
        {
            return coeffsDict_;
        }

        virtual bool sizeLocations
        (
            const pointIndexHit& hitPt,
            const vector& n,
            pointField& shapePts,
            scalarField& shapeSizes
        ) const = 0;

        //- Modify scalar argument to the cell size specified by function.
        //  Return a boolean specifying if the function was used, i.e. false if
        //  the point was not in range of the surface for a spatially varying
        //  size.
        virtual bool cellSize
        (
            const point& pt,
            scalar& size
        ) const = 0;


        virtual bool setCellSize
        (
            const pointField& pts
        )
        {
            WarningIn("cellSizeFunction::setCellSize(const pointField&)")
                << "Not overloaded."
                << endl;
            return false;
        }

        label priority() const
        {
            return priority_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
