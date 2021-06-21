/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    Foam::refinementRegions

Description
    Encapsulates queries for volume refinement ('refine all cells within
    shell').

SourceFiles
    refinementRegions.C

\*---------------------------------------------------------------------------*/

#ifndef refinementRegions_H
#define refinementRegions_H

#include "Tuple2.H"
#include "triSurfaceMesh.H"
#include "triSurfaceFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class searchableSurfaces;

/*---------------------------------------------------------------------------*\
                        Class refinementRegions Declaration
\*---------------------------------------------------------------------------*/

class refinementRegions
{
public:

    // Public data types

        //- Volume refinement controls
        enum class refineMode
        {
            inside,         // Refine all inside shell
            outside,        // ,,         outside
            distance,       // Refine based on distance to shell
            insideSpan,     // Refine based on distance internal span
            outsideSpan     // Refine based on distance external span
        };


private:

    // Private Data

        //- Reference to all geometry.
        const searchableSurfaces& allGeometry_;

        //- Indices of surfaces that are shells
        labelList shells_;

        //- Per shell whether to refine inside or outside
        List<refineMode> modes_;

        //- Per shell the list of ranges
        List<scalarField> distances_;

        //- Per shell per distance the refinement level
        labelListList levels_;

        //- Number of cells across span for refinement mode span
        labelList cellsAcrossSpan_;

        //- Surface closness field for refinement mode span
        PtrList<triSurfacePointScalarField> closeness_;


    // Private Data

        //- refineMode names
        static const NamedEnum<refineMode, 5> refineModeNames_;


    // Private Member Functions

        //- Helper function for initialisation.
        void setAndCheckLevels(const label shelli, const dictionary& dict);

        //- Specifically orient triSurfaces using a calculated point outside.
        //  Done since quite often triSurfaces not of consistent orientation
        //  which is (currently) necessary for sideness calculation
        void orient();

        //- Interpolate the surface closeness field to the given point
        scalar interpolate
        (
            const triSurfaceMesh& tris,
            const triSurfacePointScalarField& closeness,
            const point& pt,
            const label index
        ) const;

        //- Find maximum level of a shell
        void findHigherLevel
        (
            const pointField& pt,
            const label shelli,
            const scalar level0EdgeLength,
            labelList& maxLevel
        ) const;


public:

    // Constructors

        //- Construct from geometry and dictionary
        refinementRegions
        (
            const searchableSurfaces& allGeometry,
            const dictionary& shellsDict
        );


    // Member Functions

        // Query

            //- Highest shell level
            label maxLevel() const;

            //- Find shell level higher than ptLevel
            void findHigherLevel
            (
                const pointField& pt,
                const labelList& ptLevel,
                const scalar level0EdgeLength,
                labelList& maxLevel
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
