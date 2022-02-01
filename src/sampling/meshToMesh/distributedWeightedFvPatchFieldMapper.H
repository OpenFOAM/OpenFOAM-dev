/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2022 OpenFOAM Foundation
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
    // Private Data

        const label singlePatchProc_;

        const mapDistributeBase* distMapPtr_;

        const labelListList& addressing_;

        const scalarListList& weights_;

        bool hasUnmapped_;


    // Private Member Functions

        label size() const
        {
            if (singlePatchProc_ == -1)
            {
                return distMapPtr_->constructSize();
            }
            else
            {
                return addressing().size();
            }
        }

        template<class Type>
        void map(Field<Type>& f, const Field<Type>& mapF) const;

        template<class Type>
        tmp<Field<Type>> map(const Field<Type>& f) const;


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

        virtual bool hasUnmapped() const
        {
            return hasUnmapped_;
        }

        virtual const labelListList& addressing() const;

        virtual const scalarListList& weights() const;


    // Member Operators

        //- Map a field
        FOR_ALL_FIELD_TYPES(DEFINE_FIELD_MAPPER_OPERATOR, );

        //- Map a temporary field
        template<class Type>
        void operator()(Field<Type>& f, const tmp<Field<Type>>& tmapF) const;

        //- Map a temporary field
        template<class Type>
        tmp<Field<Type>> operator()(const tmp<Field<Type>>& tmapF) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "distributedWeightedFvPatchFieldMapperTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
