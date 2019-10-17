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
    Foam::Function1s::TableFile

Description
    Templated table container function where data is read from file. Data is
    stored as a list of Tuple2's. The first column is scalar. Data is read in
    the form specified by the format entry, which defaults to the native
    OpenFOAM format.

    Usage:
    \verbatim
        <entryName> tableFile;
        <entryName>Coeffs
        {
            file                dataFile;    // name of data file
            format              foam;        // data format (optional)
            outOfBounds         clamp;       // optional out-of-bounds handling
            interpolationScheme linear;      // optional interpolation method
        }
    \endverbatim

SourceFiles
    TableFile.C

See also
    FoamTableReader.C
    CsvTableReader.C

\*---------------------------------------------------------------------------*/

#ifndef TableFile_H
#define TableFile_H

#include "TableBase.H"
#include "TableReader.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1s
{

/*---------------------------------------------------------------------------*\
                           Class TableFile Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class TableFile
:
    public TableBase<Type, TableFile<Type>>
{
    // Private Data

        //- File name for table
        const fileName fName_;

        //- Table reader
        const autoPtr<TableReader<Type>> reader_;


public:

    //- Runtime type information
    TypeName("tableFile");


    // Constructors

        //- Construct from entry name and dictionary
        TableFile(const word& entryName, const dictionary& dict);

        //- Copy constructor
        TableFile(const TableFile<Type>& tbl);


    //- Destructor
    virtual ~TableFile();


    // Access

        //- Return the file name
        inline const fileName& fName() const
        {
            return fName_;
        }


    // I/O

        //- Write entries only in dictionary format
        virtual void writeEntries(Ostream& os) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const TableFile<Type>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1s
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "TableFile.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
