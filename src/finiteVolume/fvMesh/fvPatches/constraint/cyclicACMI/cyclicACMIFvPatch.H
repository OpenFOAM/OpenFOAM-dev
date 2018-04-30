/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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
    Foam::cyclicACMIFvPatch

Description
    Cyclic patch for Arbitrarily Coupled Mesh Interface (ACMI)

SourceFiles
    cyclicACMIFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicACMIFvPatch_H
#define cyclicACMIFvPatch_H

#include "cyclicAMIFvPatch.H"
#include "cyclicACMIPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class cyclicACMIFvPatch Declaration
\*---------------------------------------------------------------------------*/

class cyclicACMIFvPatch
:
    public cyclicAMIFvPatch
{
private:

    // Private member functions

        //- Update the patch areas after AMI update
        void updateAreas() const;


protected:

    // Protected Member functions

        //- Make patch weighting factors
        virtual void makeWeights(scalarField&) const;


public:

    //- Runtime type information
    TypeName(cyclicACMIPolyPatch::typeName_());


    // Constructors

        //- Construct from polyPatch
        cyclicACMIFvPatch(const polyPatch& patch, const fvBoundaryMesh& bm)
        :
           cyclicAMIFvPatch(patch, bm)
        {}


    // Member functions

        // Access

            //- Return a reference to the AMI interpolators
            virtual const PtrList<AMIPatchToPatchInterpolation>& AMIs() const
            {
                const PtrList<AMIPatchToPatchInterpolation>& AMIs =
                    cyclicAMIFvPatch::AMIs();

                updateAreas();

                return AMIs;
            }

            //- Return a reference to the AMI transformations
            virtual const List<vectorTensorTransform>& AMITransforms() const
            {
                const List<vectorTensorTransform>& AMITransforms =
                    cyclicAMIFvPatch::AMITransforms();

                updateAreas();

                return AMITransforms;
            }

            //- Return the poly patch
            const cyclicACMIPolyPatch& cyclicACMIPatch() const
            {
                return refCast<const cyclicACMIPolyPatch>(this->patch());
            }

            //- Return the neighbour patch
            const cyclicACMIFvPatch& neighbFvPatch() const
            {
                return refCast<const cyclicACMIFvPatch>
                (
                    this->boundaryMesh()[this->cyclicAMIPatch().neighbPatchID()]
                );
            }

            //- Return the non-overlap patch
            const fvPatch& nonOverlapFvPatch() const
            {
                return this->boundaryMesh()
                [
                    cyclicACMIPatch().nonOverlapPatchID()
                ];
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
