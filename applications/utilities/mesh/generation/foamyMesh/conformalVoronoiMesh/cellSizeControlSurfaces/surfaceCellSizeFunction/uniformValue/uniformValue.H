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
    Foam::uniformValue

Description

SourceFiles
    uniformValue.C

\*---------------------------------------------------------------------------*/

#ifndef uniformValue_H
#define uniformValue_H

#include "surfaceCellSizeFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class uniformValue Declaration
\*---------------------------------------------------------------------------*/

class uniformValue
:
    public surfaceCellSizeFunction
{

private:

    // Private Data

        scalar surfaceCellSize_;


public:

    //- Runtime type information
    TypeName("uniformValue");


    // Constructors

        //- Construct from components
        uniformValue
        (
            const dictionary& cellSizeFunctionDict,
            const searchableSurface& surface,
            const scalar& defaultCellSize
        );


    //- Destructor
    virtual ~uniformValue()
    {}


    // Member Functions

        virtual scalar interpolate
        (
            const point& pt,
            const label index
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
