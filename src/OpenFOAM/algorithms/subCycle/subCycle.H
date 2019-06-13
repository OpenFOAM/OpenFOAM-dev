/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::subCycle

Description
    Perform a subCycleTime on a field or list of fields.

\*---------------------------------------------------------------------------*/

#ifndef subCycle_H
#define subCycle_H

#include "subCycleTime.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class subCycleField Declaration
\*---------------------------------------------------------------------------*/

template<class GeometricField>
class subCycleField
{
    // Private Data

        //- Reference to the field being sub-cycled
        GeometricField& gf_;

        //- Reference to the field old-time field being sub-cycled
        //  Needed to avoid calls to oldTime() which may cause
        //  unexpected updates of the old-time field
        GeometricField& gf0_;

        //- Copy of the "real" old-time value of the field
        GeometricField gf_0_;


public:

    //- Type of the field
    typedef GeometricField FieldsType;

    // Constructors

        //- Construct from field and number of sub-cycles
        subCycleField(GeometricField& gf)
        :
            gf_(gf),
            gf0_(gf.oldTime()),
            gf_0_(gf0_.name() + "_", gf0_)
        {}


    //- Destructor
    ~subCycleField()
    {
        // Reset the old-time field
        gf0_ = gf_0_;

        // Correct the time index of the field to correspond to the global time
        gf_.timeIndex() = time().timeIndex();
        gf0_.timeIndex() = time().timeIndex();
    }

    //- Access to time
    const Time& time() const
    {
        return gf_.time();
    }

    //- Correct the time index of the field to correspond to
    //  the sub-cycling time.
    //
    //  The time index is incremented to protect the old-time value from
    //  being updated at the beginning of the time-loop in the case of
    //  outer iteration
    void updateTimeIndex()
    {
        gf_.timeIndex() = time().timeIndex() + 1;
        gf0_.timeIndex() = time().timeIndex() + 1;
    }
};


/*---------------------------------------------------------------------------*\
                         Class subCycleFields Declaration
\*---------------------------------------------------------------------------*/

template<class GeometricField>
class subCycleFields
{
    // Private Data

        //- List of pointers to the fields being sub-cycled
        List<GeometricField*> gfPtrs_;

        //- List of pointers to the fields old-time field being sub-cycled
        //  Needed to avoid calls to oldTime() which may cause
        //  unexpected updates of the old-time field
        List<GeometricField*> gf0Ptrs_;

        //- Copy of the "real" old-time value of the fields
        PtrList<GeometricField> gf_0Ptrs_;


public:

    //- Type of the list of fields
    typedef List<GeometricField*> FieldsType;

    // Constructors

        //- Construct from field list and number of sub-cycles
        subCycleFields(List<GeometricField*>& gfPtrs)
        :
            gfPtrs_(gfPtrs),
            gf0Ptrs_(gfPtrs.size()),
            gf_0Ptrs_(gfPtrs.size())
        {
            forAll(gfPtrs_, i)
            {
                gf0Ptrs_[i] = &gfPtrs_[i]->oldTime();

                gf_0Ptrs_.set
                (
                    i,
                    new volScalarField
                    (
                        gf0Ptrs_[i]->name() + "_",
                        *gf0Ptrs_[i]
                    )
                );
            }
        }


    //- Destructor
    ~subCycleFields()
    {
        forAll(gfPtrs_, i)
        {
            // Reset the old-time fields
            *gf0Ptrs_[i] = gf_0Ptrs_[i];

            // Correct the time index of the fields
            // to correspond to the global time
            gfPtrs_[i]->timeIndex() = time().timeIndex();
            gf0Ptrs_[i]->timeIndex() = time().timeIndex();
        }
    }

    //- Access to time
    const Time& time() const
    {
        return gfPtrs_[0]->time();
    }

    //- Correct the time index of the fields to correspond to
    //  the sub-cycling time.
    //
    //  The time index is incremented to protect the old-time value from
    //  being updated at the beginning of the time-loop in the case of
    //  outer iteration
    void updateTimeIndex()
    {
        forAll(gfPtrs_, i)
        {
            gfPtrs_[i]->timeIndex() = time().timeIndex() + 1;
            gf0Ptrs_[i]->timeIndex() = time().timeIndex() + 1;
        }
    }
};


/*---------------------------------------------------------------------------*\
                          Class subCycle Declaration
\*---------------------------------------------------------------------------*/

template
<
    class GeometricField,
    template<class> class SubCycleField = subCycleField
>
class subCycle
:
    public SubCycleField<GeometricField>,
    public subCycleTime
{

public:

    // Constructors

        //- Construct field and number of sub-cycles
        subCycle
        (
            typename SubCycleField<GeometricField>::FieldsType& gf,
            const label nSubCycles
        )
        :
            SubCycleField<GeometricField>(gf),
            subCycleTime(const_cast<Time&>(this->time()), nSubCycles)
        {
            // Update the field time index to correspond to the sub-cycle time
            this->updateTimeIndex();
        }

        //- Disallow default bitwise copy construction
        subCycle(const subCycle<GeometricField>&) = delete;


    //- Destructor
    //  End the subCycleTime, which restores the time state
    ~subCycle()
    {
        endSubCycle();
    }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const subCycle<GeometricField>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
