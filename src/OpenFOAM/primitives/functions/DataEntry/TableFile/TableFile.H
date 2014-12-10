/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::TableFile

Description
    Templated table container data entry where data is read from file.

    \verbatim
        <entryName>   tableFile;
        tableFileCoeffs
        {
            dimensions          [0 0 1 0 0]; // optional dimensions
            fileName            dataFile;    // name of data file
            outOfBounds         clamp;       // optional out-of-bounds handling
            interpolationScheme linear;      // optional interpolation method
        }
    \endverbatim

    Items are stored in a list of Tuple2's. First column is always stored as
    scalar entries.  Data is read in the form, e.g. for an entry \<entryName\>
    that is (scalar, vector):
    \verbatim
        (
            0.0 (1 2 3)
            1.0 (4 5 6)
        );
    \endverbatim


SourceFiles
    TableFile.C

\*---------------------------------------------------------------------------*/

#ifndef TableFile_H
#define TableFile_H

#include "DataEntry.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
class TableFile;

template<class Type>
Ostream& operator<<
(
    Ostream&,
    const TableFile<Type>&
);

/*---------------------------------------------------------------------------*\
                           Class TableFile Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class TableFile
:
    public DataEntry<Type>,
    public TableBase<Type>
{
    // Private data

        //- File name for csv table (optional)
        fileName fName_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const TableFile<Type>&);


public:

    //- Runtime type information
    TypeName("tableFile");


    // Constructors

        //- Construct from entry name and Istream
        TableFile(const word& entryName, const dictionary& dict);

        //- Copy constructor
        TableFile(const TableFile<Type>& tbl);

        //- Construct and return a clone
        virtual tmp<DataEntry<Type> > clone() const
        {
            return tmp<DataEntry<Type> >(new TableFile<Type>(*this));
        }


    //- Destructor
    virtual ~TableFile();


    // Member Functions

        // Manipulation

            //- Convert time
            virtual void convertTimeBase(const Time& t)
            {
                TableBase<Type>::convertTimeBase(t);
            }


        // Evaluation

            //- Return TableFile value
            virtual Type value(const scalar x) const
            {
                return TableBase<Type>::value(x);
            }

            //- Integrate between two (scalar) values
            virtual Type integrate(const scalar x1, const scalar x2) const
            {
                return TableBase<Type>::integrate(x1, x2);
            }

            //- Return dimensioned constant value
            virtual dimensioned<Type> dimValue(const scalar x) const
            {
                return TableBase<Type>::dimValue(x);
            }

            //- Integrate between two values and return dimensioned type
            virtual dimensioned<Type> dimIntegrate
            (
                const scalar x1,
                const scalar x2
            )
            {
                return TableBase<Type>::dimIntegrate(x1, x2);
            }


    // I/O

        //- Ostream Operator
        friend Ostream& operator<< <Type>
        (
            Ostream& os,
            const TableFile<Type>& tbl
        );

        //- Write in dictionary format
        virtual void writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "TableFile.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
