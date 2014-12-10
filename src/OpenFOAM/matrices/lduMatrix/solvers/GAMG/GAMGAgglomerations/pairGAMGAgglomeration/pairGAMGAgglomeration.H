/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::pairGAMGAgglomeration

Description
    Agglomerate using the pair algorithm.

SourceFiles
    pairGAMGAgglomeration.C
    pairGAMGAgglomerate.C

\*---------------------------------------------------------------------------*/

#ifndef pairGAMGAgglomeration_H
#define pairGAMGAgglomeration_H

#include "GAMGAgglomeration.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class pairGAMGAgglomeration Declaration
\*---------------------------------------------------------------------------*/

class pairGAMGAgglomeration
:
    public GAMGAgglomeration
{
    // Private data

        //- Number of levels to merge, 1 = don't merge, 2 = merge pairs etc.
        label mergeLevels_;

        //- Direction of cell loop for the current level
        static bool forward_;


protected:

    // Protected Member Functions

        //- Agglomerate all levels starting from the given face weights
        void agglomerate
        (
            const lduMesh& mesh,
            const scalarField& faceWeights
        );

        //- Disallow default bitwise copy construct
        pairGAMGAgglomeration(const pairGAMGAgglomeration&);

        //- Disallow default bitwise assignment
        void operator=(const pairGAMGAgglomeration&);


public:

    //- Runtime type information
    TypeName("pair");


    // Constructors

        //- Construct given mesh and controls
        pairGAMGAgglomeration
        (
            const lduMesh& mesh,
            const dictionary& controlDict
        );

        //- Calculate and return agglomeration
        static tmp<labelField> agglomerate
        (
            label& nCoarseCells,
            const lduAddressing& fineMatrixAddressing,
            const scalarField& faceWeights
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
