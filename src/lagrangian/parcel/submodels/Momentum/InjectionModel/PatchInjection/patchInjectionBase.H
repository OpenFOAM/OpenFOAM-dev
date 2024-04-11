/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2024 OpenFOAM Foundation
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
    Foam::patchInjectionBase

Description
    Base class for patch-based injection models.

    Class handles injecting at a random point adjacent to the patch faces to
    provide a more stochastic view of the injection process.  Patch faces are
    triangulated, and area fractions accumulated.  The fractional areas are
    then applied to determine across which face a parcel is to be injected.

SourceFiles
    PatchInjectionBase.C

\*---------------------------------------------------------------------------*/

#ifndef patchInjectionBase_H
#define patchInjectionBase_H

#include "word.H"
#include "labelList.H"
#include "scalarList.H"
#include "vectorList.H"
#include "faceList.H"
#include "triFaceList.H"
#include "barycentric.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward class declarations
class polyMesh;
class fvMesh;
class randomGenerator;

/*---------------------------------------------------------------------------*\
                     Class patchInjectionBase Declaration
\*---------------------------------------------------------------------------*/

class patchInjectionBase
{
protected:

    // Protected data

        //- Patch name
        const word patchName_;

        //- Patch ID
        const label patchId_;

        //- Cumulative sum of the total patch areas in each process
        scalarList sumProcArea_;

        //- Cumulative sum of the face areas in the patch
        scalarList sumFaceArea_;

        //- Cumulative sum of the triangle areas in each patch face
        scalarListList sumFaceTriArea_;


public:

    // Constructors

        //- Construct from mesh and patch name
        patchInjectionBase(const polyMesh& mesh, const word& patchName);

        //- Copy constructor
        patchInjectionBase(const patchInjectionBase& pib);


    //- Destructor
    virtual ~patchInjectionBase();


    // Member Functions

        //- Update patch geometry and derived info for injection locations
        virtual void topoChange(const polyMesh& mesh);

        //- Set the injection position and owner cell, tetFace and tetPt
        virtual void setPositionAndCell
        (
            const fvMesh& mesh,
            randomGenerator& rndGen,
            barycentric& coordinates,
            label& celli,
            label& tetFacei,
            label& tetPti,
            label& facei
        );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
