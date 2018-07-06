/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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
    Foam::cyclicRepeatAMIFvPatch

Description
    Repeating patch for Arbitrary Mesh Interface (AMI)

SourceFiles
    cyclicRepeatAMIFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicRepeatAMIFvPatch_H
#define cyclicRepeatAMIFvPatch_H

#include "cyclicAMIFvPatch.H"
#include "cyclicRepeatAMIPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class cyclicRepeatAMIFvPatch Declaration
\*---------------------------------------------------------------------------*/

class cyclicRepeatAMIFvPatch
:
    public cyclicAMIFvPatch
{
public:

    //- Runtime type information
    TypeName(cyclicRepeatAMIPolyPatch::typeName_());


    // Constructors

        //- Construct from polyPatch
        cyclicRepeatAMIFvPatch(const polyPatch& patch, const fvBoundaryMesh& bm)
        :
            cyclicAMIFvPatch(patch, bm)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
