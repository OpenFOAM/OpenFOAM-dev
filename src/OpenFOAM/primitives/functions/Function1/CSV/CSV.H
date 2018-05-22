/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::Function1Types::CSV

Description
    Templated CSV function.

    Reference column is always a scalar, e.g. time.

    Usage:
    \verbatim
        <entryName> csvFile;
        <entryName>Coeffs
        {
            nHeaderLine         4;          // number of header lines
            refColumn           0;          // reference column index
            componentColumns    (1 2 3);    // component column indices
            separator           ",";        // optional (defaults to ",")
            mergeSeparators     no;         // merge multiple separators
            file                "fileXYZ";  // name of csv data file
            outOfBounds         clamp;      // optional out-of-bounds handling
            interpolationScheme linear;     // optional interpolation scheme
        }
    \endverbatim

SourceFiles
    CSV.C

\*---------------------------------------------------------------------------*/

#ifndef CSV_H
#define CSV_H

#include "Function1.H"
#include "TableBase.H"
#include "Tuple2.H"
#include "labelList.H"
#include "ISstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{

/*---------------------------------------------------------------------------*\
                           Class CSV Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class CSV
:
    public TableBase<Type>
{
    // Private data

        //- Number header lines
        label nHeaderLine_;

        //- Column of the time
        label refColumn_;

        //- Labels of the components
        labelList componentColumns_;

        //- Separator character
        char separator_;

        //- Merge separators flag, e.g. ',,,' becomes ','
        bool mergeSeparators_;

        //- File name for csv table
        fileName fName_;


    // Private Member Functions

        //- Read csv data table
        void read();

        //- Read the next value from the split string
        Type readValue(const List<string>&);

        //- Disallow default bitwise assignment
        void operator=(const CSV<Type>&);


public:

    //- Runtime type information
    TypeName("csvFile");


    // Constructors

        //- Construct from entry name and dictionary
        CSV
        (
            const word& entryName,
            const dictionary& dict
        );

        //- Copy constructor
        CSV(const CSV<Type>& tbl);


    //- Destructor
    virtual ~CSV();


    // Member Functions

        //- Return const access to the file name
        virtual const fileName& fName() const;

        //- Write in dictionary format
        virtual void writeData(Ostream& os) const;
};


template<>
label CSV<label>::readValue(const List<string>& split);

template<>
Foam::scalar CSV<scalar>::readValue(const List<string>& split);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CSV.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
