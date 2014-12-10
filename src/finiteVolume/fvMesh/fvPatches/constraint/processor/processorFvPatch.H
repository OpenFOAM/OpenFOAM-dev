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
    Foam::processorFvPatch

Description
    Processor patch.

SourceFiles
    processorFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef processorFvPatch_H
#define processorFvPatch_H

#include "coupledFvPatch.H"
#include "processorLduInterface.H"
#include "processorPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class processorFvPatch Declaration
\*---------------------------------------------------------------------------*/

class processorFvPatch
:
    public coupledFvPatch,
    public processorLduInterface
{
    // Private Data

        const processorPolyPatch& procPolyPatch_;


protected:

    // Protected Member functions

        //- Make patch weighting factors
        void makeWeights(scalarField&) const;


public:

    //- Runtime type information
    TypeName(processorPolyPatch::typeName_());


    // Constructors

        //- Construct from components
        processorFvPatch(const polyPatch& patch, const fvBoundaryMesh& bm)
        :
            coupledFvPatch(patch, bm),
            procPolyPatch_(refCast<const processorPolyPatch>(patch))
        {}


    // Member functions

        //- Return communicator used for comms
        virtual int comm() const
        {
            return procPolyPatch_.comm();
        }

        //- Return processor number
        virtual int myProcNo() const
        {
            return procPolyPatch_.myProcNo();
        }

        //- Return neigbour processor number
        virtual int neighbProcNo() const
        {
            return procPolyPatch_.neighbProcNo();
        }

        //- Return message tag used for sending
        virtual int tag() const
        {
            return UPstream::msgType();
        }

        //- Return true if running parallel
        virtual bool coupled() const
        {
            if (Pstream::parRun())
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        const processorPolyPatch& procPolyPatch() const
        {
            return procPolyPatch_;
        }

        //- Are the cyclic planes parallel
        virtual bool parallel() const
        {
            return procPolyPatch_.parallel();
        }

        //- Return face transformation tensor
        virtual const tensorField& forwardT() const
        {
            return procPolyPatch_.forwardT();
        }

        //- Return neighbour-cell transformation tensor.
        virtual const tensorField& reverseT() const
        {
            return procPolyPatch_.reverseT();
        }

        //- Return delta (P to N) vectors across coupled patch
        virtual tmp<vectorField> delta() const;


        // Interface transfer functions

            //- Return the values of the given internal data adjacent to
            //  the interface as a field
            virtual tmp<labelField> interfaceInternalField
            (
                const labelUList& internalData
            ) const;

            //- Initialise neighbour field transfer
            virtual void initInternalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& internalData
            ) const;

            //- Return neighbour field
            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& internalData
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
