/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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
    Foam::cyclicRepeatAMIPointPatch

Description
    Repeat AMI point patch - place holder only

SourceFiles
    cyclicRepeatAMIPointPatch.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicRepeatAMIPointPatch_H
#define cyclicRepeatAMIPointPatch_H

#include "cyclicAMIPointPatch.H"
#include "cyclicRepeatAMIPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class cyclicRepeatAMIPointPatch Declaration
\*---------------------------------------------------------------------------*/

class cyclicRepeatAMIPointPatch
:
    public cyclicAMIPointPatch
{
public:

    //- Runtime type information
    TypeName(cyclicRepeatAMIPolyPatch::typeName_());


    // Constructors

        //- Construct from components
        cyclicRepeatAMIPointPatch
        (
            const polyPatch& patch,
            const pointBoundaryMesh& bm
        );


    //- Destructor
    virtual ~cyclicRepeatAMIPointPatch();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
