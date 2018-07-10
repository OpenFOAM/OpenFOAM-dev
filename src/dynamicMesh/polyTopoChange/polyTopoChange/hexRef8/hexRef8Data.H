/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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
    Foam::hexRef8Data

Description
    Various for reading/decomposing/reconstructing/distributing refinement
    data.

SourceFiles
    hexRef8Data.C

\*---------------------------------------------------------------------------*/

#ifndef hexRef8Data_H
#define hexRef8Data_H

#include "labelIOList.H"
#include "uniformDimensionedFields.H"
#include "UPtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class mapPolyMesh;
class mapDistributePolyMesh;
class refinementHistory;
class fvMesh;

/*---------------------------------------------------------------------------*\
                           Class hexRef8Data Declaration
\*---------------------------------------------------------------------------*/

class hexRef8Data
{

private:

    // Private data

        autoPtr<labelIOList> cellLevelPtr_;

        autoPtr<labelIOList> pointLevelPtr_;

        autoPtr<uniformDimensionedScalarField> level0EdgePtr_;

        autoPtr<refinementHistory> refHistoryPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        hexRef8Data(const hexRef8Data&);

        //- Disallow default bitwise assignment
        void operator=(const hexRef8Data&);


public:

    // Constructors

        //- Construct read. Has special provision for only some processors
        //  having the files so can be used in redistribution.
        hexRef8Data(const IOobject& io);

        //- Construct as subset
        hexRef8Data
        (
            const IOobject& io,
            const hexRef8Data&,
            const labelList& cellMap,
            const labelList& pointMap
        );

        //- Construct from multiple hexRef8Data
        hexRef8Data
        (
            const IOobject& io,
            const UPtrList<const labelList>& cellMaps,
            const UPtrList<const labelList>& pointMaps,
            const UPtrList<const hexRef8Data>&
        );


    //- Destructor
    ~hexRef8Data();


    // Member Functions

        //- Parallel synchronise. This enforces valid objects on all processors
        //  (even if they don't have a mesh). Used by redistributePar.
        void sync(const IOobject& io);

        //- In-place update for topology changes
        void updateMesh(const mapPolyMesh&);

        //- In-place distribute
        void distribute(const mapDistributePolyMesh&);

        //- Write
        bool write() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
