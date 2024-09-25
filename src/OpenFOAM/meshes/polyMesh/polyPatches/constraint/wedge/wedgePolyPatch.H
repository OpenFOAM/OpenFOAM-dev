/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
    Foam::wedgePolyPatch

Description
    Wedge front and back plane patch.

SourceFiles
    wedgePolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef wedgePolyPatch_H
#define wedgePolyPatch_H

#include "polyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class wedgePolyPatch Declaration
\*---------------------------------------------------------------------------*/

class wedgePolyPatch
:
    public polyPatch
{
    // Private Data

        //- Axis of the wedge
        vector axis_;

        //- Centre normal between the wedge boundaries
        vector centreNormal_;

        //- Normal to the patch
        vector n_;

        //- Cosine of the wedge angle
        scalar cosAngle_;

        //- Face transformation tensor
        tensor faceT_;

        //- Neighbour-cell transformation tensor
        tensor cellT_;

        //- Index of the opposite wedge patch
        mutable label oppositePatchIndex_;


    // Private Member Functions

        //- Return whether the given face is a (non-degenerate) quad
        bool isQuadFace(const label) const;

        //- Return whether the given face is part of a wedge patch
        bool isWedgeFace(const label) const;

        //- Return the wedge face opposite the given wedge face
        label oppositeWedgeFace(const label) const;


protected:

    // Protected Member Functions

        //- Calculate the patch geometry
        virtual void calcGeometry(PstreamBuffers&);


public:

    //- Runtime type information
    TypeName("wedge");


    // Constructors

        //- Construct from components
        wedgePolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType
        );

        //- Construct from dictionary
        wedgePolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType
        );

        //- Construct as copy, resetting the boundary mesh
        wedgePolyPatch(const wedgePolyPatch&, const polyBoundaryMesh&);

        //- Construct given the original patch and resetting the
        //  face list and boundary mesh information
        wedgePolyPatch
        (
            const wedgePolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        );

        //- Construct and return a clone, resetting the boundary mesh
        virtual autoPtr<polyPatch> clone(const polyBoundaryMesh& bm) const
        {
            return autoPtr<polyPatch>(new wedgePolyPatch(*this, bm));
        }

        //- Construct and return a clone, resetting the face list
        //  and boundary mesh
        virtual autoPtr<polyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        ) const
        {
            return autoPtr<polyPatch>
            (
                new wedgePolyPatch(*this, bm, index, newSize, newStart)
            );
        }


    // Member Functions

        // Access

            //- Return axis of the wedge
            const vector& axis() const
            {
                return axis_;
            }

            //- Return plane normal between the wedge boundaries
            const vector& centreNormal() const
            {
                return centreNormal_;
            }

            //- Return the normal to the patch
            const vector& n() const
            {
                return n_;
            }

            //- Return the cosine of the wedge angle
            scalar cosAngle() const
            {
                return cosAngle_;
            }

            //- Return face transformation tensor
            const tensor& faceT() const
            {
                return faceT_;
            }

            //- Return neighbour-cell transformation tensor
            const tensor& cellT() const
            {
                return cellT_;
            }

            //- Return the index of the opposite wedge patch
            label oppositePatchIndex() const;

            //- Return a reference to the opposite wedge patch
            const wedgePolyPatch& oppositePatch() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
