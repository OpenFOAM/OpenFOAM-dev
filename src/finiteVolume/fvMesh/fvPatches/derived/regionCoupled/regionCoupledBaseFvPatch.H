/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::regionCoupledBaseFvPatch

Description
    Base class of regionCoupledFvPatch with common functionality for
    regionCoupledFvPatch and regionCoupledWallFvPatch

SourceFiles
    regionCoupledBaseFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef regionCoupledBaseFvPatch_H
#define regionCoupledBaseFvPatch_H

#include "regionCoupledBase.H"
#include "fvMesh.H"
#include "Time.H"
#include "polyPatch.H"
#include "regionCoupledLduInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class regionCoupledBaseFvPatch Declaration
\*---------------------------------------------------------------------------*/

class regionCoupledBaseFvPatch
:
    public regionCoupledLduInterface
{
protected:

    // Protected data

        //- Reference to regionCoupledBase
        const regionCoupledBase& regionCoupledBase_;

        //- Reference to fvPatch
        const fvPatch& patch_;


    // Protected members

        //- Returns fvMesh
        const fvMesh& nbrFvMesh() const
        {
            return
            (
                patch_.boundaryMesh().mesh().time().lookupObject<fvMesh>
                (
                    regionCoupledBase_.nbrRegionName()
                )
            );
        }


public:

    //- Runtime type information
    TypeName("regionCoupledBase");


    // Constructors

        //- Construct from polyPatch
        regionCoupledBaseFvPatch
        (
            const polyPatch& pp,
            const fvPatch& patch
        )
        :
            regionCoupledLduInterface(),
            regionCoupledBase_
            (
                refCast<const regionCoupledBase>(pp)
            ),
            patch_(patch)
        {}


    //- Destructor
    virtual ~regionCoupledBaseFvPatch()
    {}


    // Member Functions

        // Access

            //- Return neighbour
            virtual label neighbPatchID() const
            {
                return regionCoupledBase_.neighbPatchID();
            }

            //- Is it the owner?
            virtual bool owner() const
            {
                return regionCoupledBase_.owner();
            }

            //- Return regionCoupledBase neighb Patch
            virtual const regionCoupledBaseFvPatch& neighbPatch() const
            {
                return refCast<const regionCoupledBaseFvPatch>
                (
                    nbrFvMesh().boundary()
                    [
                        regionCoupledBase_.neighbPatchID()
                    ]
                );
            }

            //- Return a reference to the AMI interpolator
            virtual const AMIPatchToPatchInterpolation& AMI() const
            {
                return regionCoupledBase_.AMI();
            }

            //- Returns neighbour polyMesh
            virtual const polyMesh& nbrMesh() const
            {
                return
                (
                    patch_.boundaryMesh().mesh().time().lookupObject<polyMesh>
                    (
                        regionCoupledBase_.nbrRegionName()
                    )
                );
            }

            //- Return fvPatch
            const fvPatch& patch() const
            {
                return patch_;
            }

            //- Returns if it is the same Region
            bool sameRegion() const
            {
                return regionCoupledBase_.sameRegion();
            }

            //- Return regionCoupledPolyPatch
            const regionCoupledBase& regionCoupledPatch() const
            {
                return regionCoupledBase_;
            }

            //- Return neighbor fvPatch
            const fvPatch& neighbFvPatch() const
            {
                return refCast<const fvPatch>
                (
                    nbrFvMesh().boundary()
                    [
                        regionCoupledBase_.neighbPatchID()
                    ]
                );
            }

            //- Return the interface type
            const word& regionCoupleType() const
            {
                return regionCoupledBase_.regionCoupleType();
            }


             //- Return faceCell addressing
            virtual const labelUList& faceCells() const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
