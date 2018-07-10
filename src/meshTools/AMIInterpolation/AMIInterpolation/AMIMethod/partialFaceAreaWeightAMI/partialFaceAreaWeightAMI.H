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
    Foam::partialFaceAreaWeightAMI

Description
    Partial face area weighted Arbitrary Mesh Interface (AMI) method

SourceFiles
    partialFaceAreaWeightAMI.C

\*---------------------------------------------------------------------------*/

#ifndef partialFaceAreaWeightAMI_H
#define partialFaceAreaWeightAMI_H

#include "faceAreaWeightAMI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class partialFaceAreaWeightAMI Declaration
\*---------------------------------------------------------------------------*/

template<class SourcePatch, class TargetPatch>
class partialFaceAreaWeightAMI
:
    public faceAreaWeightAMI<SourcePatch, TargetPatch>
{

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        partialFaceAreaWeightAMI(const partialFaceAreaWeightAMI&);

        //- Disallow default bitwise assignment
        void operator=(const partialFaceAreaWeightAMI&);

        // Marching front

            //- Set the source and target seed faces
            virtual void setNextFaces
            (
                label& startSeedI,
                label& srcFacei,
                label& tgtFacei,
                const boolList& mapFlag,
                labelList& seedFaces,
                const DynamicList<label>& visitedFaces,
                bool errorOnNotFound = true
            ) const;


public:

    //- Runtime type information
    TypeName("partialFaceAreaWeightAMI");


    // Constructors

        //- Construct from components
        partialFaceAreaWeightAMI
        (
            const SourcePatch& srcPatch,
            const TargetPatch& tgtPatch,
            const scalarField& srcMagSf,
            const scalarField& tgtMagSf,
            const faceAreaIntersect::triangulationMode& triMode,
            const bool reverseTarget = false,
            const bool requireMatch = true
        );


    //- Destructor
    virtual ~partialFaceAreaWeightAMI();


    // Member Functions

        // Access

            //- Flag to indicate that interpolation patches are conformal
            virtual bool conformal() const;


        // Manipulation

            //- Update addressing and weights
            virtual void calculate
            (
                labelListList& srcAddress,
                scalarListList& srcWeights,
                labelListList& tgtAddress,
                scalarListList& tgtWeights,
                label srcFacei = -1,
                label tgtFacei = -1
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "partialFaceAreaWeightAMI.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
