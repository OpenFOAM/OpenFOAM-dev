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
    Foam::distributedWeightedFvPatchFieldMapper

Description
    FieldMapper with weighted mapping from (optionally remote) quantities.

\*---------------------------------------------------------------------------*/

#ifndef distributedWeightedFvPatchFieldMapper_H
#define distributedWeightedFvPatchFieldMapper_H

#include "fvPatchFieldMapper.H"
#include "mapDistributeBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class distributedWeightedFvPatchFieldMapper Declaration
\*---------------------------------------------------------------------------*/

class distributedWeightedFvPatchFieldMapper
:
    public fvPatchFieldMapper
{
    const label singlePatchProc_;

    const mapDistributeBase* distMapPtr_;

    const labelListList& addressing_;

    const scalarListList& weights_;

    bool hasUnmapped_;

public:

    // Constructors

        //- Construct given addressing
        distributedWeightedFvPatchFieldMapper
        (
            const label singlePatchProc,
            const mapDistributeBase* distMapPtr,
            const labelListList& addressing,
            const scalarListList& weights
        )
        :
            singlePatchProc_(singlePatchProc),
            distMapPtr_(distMapPtr),
            addressing_(addressing),
            weights_(weights),
            hasUnmapped_(false)
        {
            forAll(addressing_, i)
            {
                if (addressing_[i].size() == 0)
                {
                    hasUnmapped_ = true;
                }
            }

            if ((singlePatchProc_ == -1) != (distMapPtr_ != nullptr))
            {
                FatalErrorIn
                (
                    "distributedWeightedFvPatchFieldMapper::"
                    "distributedWeightedFvPatchFieldMapper(..)"
                )   << "Supply a mapDistributeBase if and only if "
                    << "singlePatchProc is -1"
                    << " singlePatchProc_:" << singlePatchProc_
                    << " distMapPtr_:" << (distMapPtr_ != nullptr)
                    << exit(FatalError);
            }
        }

    //- Destructor
    virtual ~distributedWeightedFvPatchFieldMapper()
    {}


    // Member Functions

        virtual label size() const
        {
            if (distributed())
            {
                return distributeMap().constructSize();
            }
            else
            {
                return addressing().size();
            }
        }

        virtual bool direct() const
        {
            return false;
        }

        virtual bool distributed() const
        {
            return singlePatchProc_ == -1;
        }

        virtual const mapDistributeBase& distributeMap() const
        {
            if (!distMapPtr_)
            {
                FatalErrorIn
                (
                    "distributedWeightedFvPatchFieldMapper::"
                    "distributeMap()"
                )   << "Cannot ask for distributeMap on a non-distributed"
                    << " mapper" << exit(FatalError);
            }
            return *distMapPtr_;
        }

        virtual bool hasUnmapped() const
        {
            return hasUnmapped_;
        }

        virtual const labelListList& addressing() const
        {
            return addressing_;
        }

        virtual const scalarListList& weights() const
        {
            return weights_;
        }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
