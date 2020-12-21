/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::Function1s::Table

Description
    Templated interpolated tabulated data Function1.

    Items are stored in a list of Tuple2's. First column is always stored as
    scalar entries. Data is read in Tuple2 form:

    Usage:
    \verbatim
        <name>   table
        (
            (0.0 (1 2 3))
            (1.0 (4 5 6))
        );
    \endverbatim

    or in dictionary form which supports the setting of options, e.g.
    \verbatim
        <name>   table;

        values
        (
            (0.0 (1 2 3))
            (1.0 (4 5 6))
        );

        outOfBounds         clamp;     // optional out-of-bounds handling
        interpolationScheme linear;    // optional interpolation method
    \endverbatim

    or in sub-dictionary form which avoids clashes between table entries and
    other entries in the dictionary:

    \verbatim
        <name>
        {
            type    table;

            values
            (
                (0.0 (1 2 3))
                (1.0 (4 5 6))
            );

            outOfBounds         clamp;     // optional out-of-bounds handling
            interpolationScheme linear;    // optional interpolation method
        }
    \endverbatim

    The data may be read from a separate file in either native or CSV format:

    Usage:
    \verbatim
        <name>
        {
            type                table;
            file                "<file path>"; // Name/path of thedata file
            format              foam;      // data format (optional)
            outOfBounds         clamp;     // optional out-of-bounds handling
            interpolationScheme linear;    // optional interpolation method
        }
    \endverbatim

SourceFiles
    Table.C

See also
    FoamTableReader.C
    CsvTableReader.C

\*---------------------------------------------------------------------------*/

#ifndef Table_H
#define Table_H

#include "tableBase.H"
#include "Function1.H"
#include "TableReader.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class interpolationWeights;

namespace Function1s
{

/*---------------------------------------------------------------------------*\
                          Class Table Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Table
:
    public tableBase,
    public FieldFunction1<Type, Table<Type>>
{
    // Private Data

        //- Enumeration for handling out-of-bound values
        const tableBase::boundsHandling boundsHandling_;

        //- Interpolation type
        const word interpolationScheme_;

        //- Table data
        List<Tuple2<scalar, Type>> values_;

        //- Extracted values
        mutable autoPtr<scalarField> tableSamplesPtr_;

        //- Interpolator method
        mutable autoPtr<interpolationWeights> interpolatorPtr_;

        //- Cached indices
        mutable labelList indices_;

        //- Cached weights
        mutable scalarField weights_;

        //- Table reader
        const autoPtr<TableReader<Type>> reader_;


    // Private Member Functions

        //- Return (demand driven) interpolator
        const interpolationWeights& interpolator() const;

        //- Check the table for size and consistency
        void check() const;

        //- Bound the argument to the table. Errors or warns, or shifts the
        //  value if the table repeats. Does not clamp to the ends of the table
        //  as the interpolator already performs that function.
        scalar bound(const scalar x) const;


public:

    //- Runtime type information
    TypeName("table");


    // Constructors

        //- Construct from components
        Table
        (
            const word& name,
            const tableBase::boundsHandling boundsHandling,
            const word& interpolationScheme,
            const List<Tuple2<scalar, Type>>& table
        );

        //- Construct from name and dictionary
        Table(const word& name, const dictionary& dict);

        //- Copy constructor
        Table(const Table<Type>& tbl);


    //- Destructor
    virtual ~Table();


    // Member Functions

        //- Return Table value as a function of scalar x
        virtual Type value(const scalar x) const;

        //- Integrate between two scalars
        virtual Type integral(const scalar x1, const scalar x2) const;

        //- Return the reference values
        virtual tmp<scalarField> x() const;

        //- Return the dependent values
        virtual tmp<Field<Type>> y() const;

        //- Write data to dictionary stream
        virtual void write(Ostream& os) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const Table<Type>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1s
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Table.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
