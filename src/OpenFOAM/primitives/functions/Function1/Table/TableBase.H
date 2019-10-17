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
    Foam::Function1Types::TableBase

Description
    Base class for table with bounds handling, interpolation and integration

SourceFiles
    TableBase.C

\*---------------------------------------------------------------------------*/

#ifndef TableBase_H
#define TableBase_H

#include "tableBase.H"
#include "Function1.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class interpolationWeights;

namespace Function1Types
{

/*---------------------------------------------------------------------------*\
                        Class TableBase Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class TableBase
:
    public tableBase,
    public Function1<Type>
{
protected:

    // Protected data

        //- Table name
        const word name_;

        //- Enumeration for handling out-of-bound values
        const tableBase::boundsHandling boundsHandling_;

        //- Interpolation type
        const word interpolationScheme_;

        //- Table data
        List<Tuple2<scalar, Type>> table_;

        //- Extracted values
        mutable autoPtr<scalarField> tableSamplesPtr_;

        //- Interpolator method
        mutable autoPtr<interpolationWeights> interpolatorPtr_;

        //- Cached indices
        mutable labelList indices_;

        //- Cached weights
        mutable scalarField weights_;


    // Protected Member Functions

        //- Return (demand driven) interpolator
        const interpolationWeights& interpolator() const;

        //- Check the table for size and consistency
        void check() const;

        //- Bound the argument to the table. Errors or warns, or shifts the
        //  value if the table repeats. Does not clamp to the ends of the table
        //  as the interpolator already performs that function.
        scalar bound(const scalar x) const;


private:

        //- Disallow default bitwise assignment
        void operator=(const TableBase<Type>&) = delete;


public:

    // Constructors

        //- Construct from dictionary - note table is not populated
        TableBase(const word& name, const dictionary& dict);

        //- Copy constructor. Note: steals interpolator, tableSamples
        TableBase(const TableBase<Type>& tbl);


    //- Destructor
    virtual ~TableBase();


    // Member Functions

        //- Return Table value
        virtual Type value(const scalar x) const;

        //- Integrate between two (scalar) values
        virtual Type integrate(const scalar x1, const scalar x2) const;

        //- Return the reference values
        virtual tmp<scalarField> x() const;

        //- Return the dependent values
        virtual tmp<Field<Type>> y() const;

        //- Write all table data in dictionary format
        virtual void writeData(Ostream& os) const;

        //- Write keywords only in dictionary format. Used for non-inline
        //  table types
        virtual void writeEntries(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "TableBase.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
