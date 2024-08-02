/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2024 OpenFOAM Foundation
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
    Foam::mappedExtrudedPatchBase

Description
    Engine which provides mapping between two patches
    and which offsets the geometry to account for the extrusion thickness

SourceFiles
    mappedExtrudedPatchBase.C

\*---------------------------------------------------------------------------*/

#ifndef mappedExtrudedPatchBase_H
#define mappedExtrudedPatchBase_H

#include "mappedPatchBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class mappedExtrudedPatchBase Declaration
\*---------------------------------------------------------------------------*/

class mappedExtrudedPatchBase
:
    public mappedPatchBase
{
    // Private Member Data

        //- Is this the in the extruded region?
        mutable bool isExtrudedRegion_;

        //- The face areas on the bottom patch
        mutable autoPtr<vectorField> bottomFaceAreasPtr_;

        //- The face centres on the bottom patch
        mutable autoPtr<pointField> bottomFaceCentresPtr_;

        //- The local points on the bottom patch
        mutable autoPtr<pointField> bottomLocalPointsPtr_;


protected:

    // Protected Member Functions

        //- Get the face-areas for the patch
        virtual tmp<vectorField> patchFaceAreas() const;

        //- Get the face-centres for the patch
        virtual tmp<pointField> patchFaceCentres() const;

        //- Get the local points for the patch
        virtual tmp<pointField> patchLocalPoints() const;

        //- Clear out data on mesh change
        virtual void clearOut() const;


public:

    //- Runtime type information
    TypeName("mappedExtrudedPatchBase");


    // Constructors

        //- Construct from patch
        mappedExtrudedPatchBase(const polyPatch&);

        //- Construct from components
        mappedExtrudedPatchBase
        (
            const polyPatch& pp,
            const word& nbrRegionName,
            const word& nbrPatchName,
            const bool isExtrudedRegion,
            const cyclicTransform& transform
        );

        //- Construct from dictionary
        mappedExtrudedPatchBase
        (
            const polyPatch& pp,
            const dictionary& dict,
            const bool defaultTransformIsNone
        );

        //- Construct as copy, resetting the boundary mesh
        mappedExtrudedPatchBase
        (
            const polyPatch&,
            const mappedExtrudedPatchBase&
        );


    //- Destructor
    virtual ~mappedExtrudedPatchBase();


    // Member Functions

        //- Write the polyPatch data as a dictionary
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
