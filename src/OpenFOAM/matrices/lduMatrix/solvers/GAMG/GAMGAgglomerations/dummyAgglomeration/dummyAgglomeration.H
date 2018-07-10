/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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
    Foam::dummyAgglomeration

Description
    Agglomerate without combining cells. Used for testing.

SourceFiles
    dummyAgglomeration.C

\*---------------------------------------------------------------------------*/

#ifndef dummyAgglomeration_H
#define dummyAgglomeration_H

#include "GAMGAgglomeration.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dummyAgglomeration Declaration
\*---------------------------------------------------------------------------*/

class dummyAgglomeration
:
    public GAMGAgglomeration
{
    // Private data

        //- Preset number of levels
        label nLevels_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        dummyAgglomeration(const dummyAgglomeration&);

        //- Disallow default bitwise assignment
        void operator=(const dummyAgglomeration&);


public:

    //- Runtime type information
    TypeName("dummy");


    // Constructors

        //- Construct given mesh and controls
        dummyAgglomeration
        (
            const lduMesh& mesh,
            const dictionary& controlDict
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
