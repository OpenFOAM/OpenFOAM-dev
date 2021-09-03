/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2021 OpenFOAM Foundation
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
    Foam::directMethod

Description
    Direct (one-to-one cell correspondence) mesh-to-mesh interpolation class

    Volume conservative.

SourceFiles
    directMethod.C

\*---------------------------------------------------------------------------*/

#ifndef directMethod_H
#define directMethod_H

#include "meshToMeshMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class directMethod Declaration
\*---------------------------------------------------------------------------*/

class directMethod
:
    public meshToMeshMethod
{
protected:

    // Protected Member Functions

        //- Return the true if cells intersect
        virtual bool intersect
        (
            const label srcCelli,
            const label tgtCelli
        ) const;

        //- Find indices of overlapping cells in src and tgt meshes - returns
        //  true if found a matching pair
        virtual bool findInitialSeeds
        (
            const labelList& srcCellIDs,
            const boolList& mapFlag,
            const label startSeedI,
            label& srcSeedI,
            label& tgtSeedI
        ) const;

        //- Calculate the mesh-to-mesh addressing and weights
        virtual void calculateAddressing
        (
            labelListList& srcToTgtCellAddr,
            scalarListList& srcToTgtCellWght,
            labelListList& tgtToSrcCellAddr,
            scalarListList& tgtToSrcCellWght,
            const label srcSeedI,
            const label tgtSeedI,
            const labelList& srcCellIDs,
            boolList& mapFlag,
            label& startSeedI
        );

        //- Append to list of src mesh seed indices
        virtual void appendToDirectSeeds
        (
            boolList& mapFlag,
            labelList& srcTgtSeed,
            DynamicList<label>& srcSeeds,
            label& srcSeedI,
            label& tgtSeedI
        ) const;


public:

    //- Run-time type information
    TypeName("direct");

    // Constructors

        //- Construct from source and target meshes
        directMethod(const polyMesh& src, const polyMesh& tgt);

        //- Disallow default bitwise copy construction
        directMethod(const directMethod&) = delete;


    //- Destructor
    virtual ~directMethod();


    // Member Functions

        //- Return the corresponding AMI method for patch interpolation
        virtual const word& AMImethod() const;

        // Evaluate

            //- Calculate addressing and weights
            virtual void calculate
            (
                labelListList& srcToTgtAddr,
                scalarListList& srcToTgtWght,
                labelListList& tgtToTgtAddr,
                scalarListList& tgtToTgtWght
            );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const directMethod&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
