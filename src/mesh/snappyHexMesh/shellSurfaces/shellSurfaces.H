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
    Foam::shellSurfaces

Description
    Encapsulates queries for volume refinement ('refine all cells within
    shell').

SourceFiles
    shellSurfaces.C

\*---------------------------------------------------------------------------*/

#ifndef shellSurfaces_H
#define shellSurfaces_H

#include "searchableSurface.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class searchableSurfaces;

/*---------------------------------------------------------------------------*\
                           Class shellSurfaces Declaration
\*---------------------------------------------------------------------------*/

class shellSurfaces
{
public:

    // Public data types

        //- Volume refinement controls
        enum refineMode
        {
            INSIDE,         // Refine all inside shell
            OUTSIDE,        // ,,         outside
            DISTANCE        // Refine based on distance to shell
        };


private:

    // Private data

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


    // Private data

        //- refineMode names
        static const NamedEnum<refineMode, 3> refineModeNames_;


    // Private Member Functions

        //- Helper function for initialisation.
        void setAndCheckLevels
        (
            const label shellI,
            const List<Tuple2<scalar, label>>&
        );

        void orient();

        void findHigherLevel
        (
            const pointField& pt,
            const label shellI,
            labelList& maxLevel
        ) const;

public:

    // Constructors

        //- Construct from geometry and dictionary
        shellSurfaces
        (
            const searchableSurfaces& allGeometry,
            const dictionary& shellsDict
        );


    // Member Functions

        // Access

            // const List<scalarField>& distances() const
            //{
            //    return distances_;
            //}
            //
            ////- Per shell per distance the refinement level
            // const labelListList& levels() const
            //{
            //    return levels_;
            //}


        // Query

            //- Highest shell level
            label maxLevel() const;

            //- Find shell level higher than ptLevel
            void findHigherLevel
            (
                const pointField& pt,
                const labelList& ptLevel,
                labelList& maxLevel
            ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
