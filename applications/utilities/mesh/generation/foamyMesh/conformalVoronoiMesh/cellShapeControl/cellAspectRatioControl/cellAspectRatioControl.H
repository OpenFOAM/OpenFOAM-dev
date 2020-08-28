/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2020 OpenFOAM Foundation
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
    Foam::cellAspectRatioControl

Description

SourceFiles
    cellAspectRatioControl.C

\*---------------------------------------------------------------------------*/

#ifndef cellAspectRatioControl_H
#define cellAspectRatioControl_H

#include "dictionary.H"
#include "vector.H"
#include "scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                   Class cellAspectRatioControl Declaration
\*---------------------------------------------------------------------------*/

class cellAspectRatioControl
{
    // Private Data

        const dictionary aspectRatioDict_;

        const scalar aspectRatio_;

        vector aspectRatioDirection_;


public:

    // Constructors

        //- Construct from dictionary
        cellAspectRatioControl(const dictionary& motionDict);

        //- Disallow default bitwise copy construction
        cellAspectRatioControl(const cellAspectRatioControl&) = delete;


    //- Destructor
    virtual ~cellAspectRatioControl();


    // Member Functions

        // Query

            void updateCellSizeAndFaceArea
            (
                vector& alignmentDir,
                scalar& targetFaceArea,
                scalar& targetCellSize
            ) const;

            void updateDeltaVector
            (
                const vector& alignmentDir,
                const scalar targetCellSize,
                const scalar rABMag,
                vector& delta
            ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const cellAspectRatioControl&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
