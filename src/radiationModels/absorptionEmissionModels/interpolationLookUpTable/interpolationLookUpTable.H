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
    Foam::interpolationLookUpTable

Description
    A list of lists. Interpolates based on the first dimension.
    The values must be positive and monotonically increasing in each dimension

Note
    - Accessing an empty list results in an error.
    - Accessing a list with a single element always returns the same value.

SourceFiles
    interpolationLookUpTable.C

\*---------------------------------------------------------------------------*/

#ifndef interpolationLookUpTable_H
#define interpolationLookUpTable_H

#include "List.H"
#include "ListOps.H"
#include "scalarField.H"
#include "HashTable.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
namespace absorptionEmissionModels
{

/*---------------------------------------------------------------------------*\
                  Class interpolationLookUpTable Declaration
\*---------------------------------------------------------------------------*/

class interpolationLookUpTable
:
    public List<scalarField>
{
    // Private Data

        //- File name
        fileName fileName_;

        //- Table dimensions
        List<label> dim_;

        //- Min on each dimension
        List<scalar> min_;

        //- Deltas on each dimension
        List<scalar> delta_;

        //- Maximum on each dimension
        List<scalar> max_;

        //- Dictionary entries
        List<dictionary> entries_;

        //- Output dictionaries
        List<dictionary> output_;

        //- Input indices from the lookup table
        List<label> entryIndices_;

        //- Output indices from the lookup Table
        List<label> outputIndices_;

        //- Field names and indices
        HashTable<label> fieldIndices_;

        //- Output list containing input and interpolation values of outputs
        List<scalar> interpOutput_;


    // Private Member Functions

        //- Read the table of data from file
        void readTable(const word& instance, const objectRegistry&);

        //- Dimension table from dictionaries input and output
        void dimensionTable();

        //- Find table index by scalarList
        label index(const List<scalar>&, const bool lastDim=true) const;

        //- Find table index by scalar
        label index(const scalar) const;

        //- Check range of lookup value
        bool checkRange(const scalar, const label) const;

        //- Interpolate function returning a scalar
        scalar interpolate
        (
            const label lo,
            const label hi,
            const scalar lookUpValue,
            const label ofield,
            const label interfield
        ) const;

        // Check list is monotonically increasing
        void check() const;

        // find hi index, interpolate and populate interpOutput_
        void findHi(const label lo, const scalar retvals);


public:

    // Constructors

        //- Construct null
        interpolationLookUpTable();

        //- Construct given the name of the file containing the table of data
        interpolationLookUpTable
        (
            const fileName&,
            const word& instance,
            const objectRegistry&
        );

         //- Construct from dictionary
        interpolationLookUpTable(const dictionary&);

        //- Construct copy
        interpolationLookUpTable(const interpolationLookUpTable&);


    // Member Functions

        //- Return true if the field exists in the table
        bool found(const word& fieldName) const;

        //- Return the output list given a single input scalar
        const List<scalar>& lookUp(const scalar);

        //- Write lookup table to filename.
        void write
        (
            Ostream&,
            const fileName&,
            const word& instance,
            const objectRegistry&
        ) const;


    // Access

        //- Return the index of a field by name
        inline label findFieldIndex(const word& fieldName) const;

        //- Return const access to the output dictionaries
        inline const List<dictionary>& output() const;

        //- Return const access tp the dictionary entries
        inline const List<dictionary>& entries() const;

        //- Return const access to the list of min dimensions
        inline const List<scalar>& min() const;

        //- Return const access to the list of dimensions
        inline const List<label>& dim() const;

        //- Return const access to the deltas in each dimension
        inline const List<scalar>& delta() const;

        //- Return const access to the list of max dimensions
        inline const List<scalar>& max() const;

        //- Return const access to the table name
        inline word tableName() const;


     // Member Operators

        //- Return an element
        const scalarField& operator[](const label) const;

        //- Return an element
        scalarField& operator[](const label);

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace absorptionEmissionModels
} // End namespace radiationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "interpolationLookUpTableI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
