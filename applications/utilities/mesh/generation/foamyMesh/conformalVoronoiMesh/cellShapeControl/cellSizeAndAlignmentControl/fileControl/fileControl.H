/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2019 OpenFOAM Foundation
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
    Foam::fileControl

Description

SourceFiles
    fileControl.C

\*---------------------------------------------------------------------------*/

#ifndef fileControl_H
#define fileControl_H

#include "cellSizeAndAlignmentControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                      Class fileControl Declaration
\*---------------------------------------------------------------------------*/

class fileControl
:
    public cellSizeAndAlignmentControl
{
    // Private Data

        const fileName pointsFile_;

        const fileName sizesFile_;

        const fileName alignmentsFile_;

        label maxPriority_;


public:

    //- Runtime type information
    TypeName("fileControl");


    // Constructors

        //- Construct from dictionary and references to conformalVoronoiMesh and
        //  searchableSurfaces
        fileControl
        (
            const Time& runTime,
            const word& name,
            const dictionary& controlFunctionDict,
            const conformationSurfaces& geometryToConformTo,
            const scalar& defaultCellSize
        );

        //- Disallow default bitwise copy construction
        fileControl(const fileControl&) = delete;


    //- Destructor
    ~fileControl();


    // Member Functions

        // Access

            virtual label maxPriority() const
            {
                return maxPriority_;
            }


        // Edit

            virtual void cellSizeFunctionVertices
            (
                DynamicList<Foam::point>& pts,
                DynamicList<scalar>& sizes
            ) const;

            virtual void initialVertices
            (
                pointField& pts,
                scalarField& sizes,
                triadField& alignments
            ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fileControl&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
