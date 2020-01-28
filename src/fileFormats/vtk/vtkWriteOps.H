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
    Foam::vtkWriteOps

Description
    VTK ASCII and binary write functions

SourceFiles
    vtkWriteOps.C

\*---------------------------------------------------------------------------*/

#ifndef vtkWriteOps_H
#define vtkWriteOps_H

#include "floatScalar.H"
#include "point.H"
#include "tensor.H"
#include "labelList.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Namespace vtkWriteOps Declaration
\*---------------------------------------------------------------------------*/

namespace vtkWriteOps
{
    //- Swap halves of word
    void swapWord(label& word32);

    //- Swap halves of word
    void swapWords(const label nWords, label* words32);

    //- Write header
    void writeHeader
    (
        std::ostream&,
        const bool isBinary,
        const std::string& title
    );

    void writeCellDataHeader
    (
        std::ostream&,
        const label nCells,
        const label nFields
    );

    void writePointDataHeader
    (
        std::ostream&,
        const label nPoints,
        const label nFields
    );


    //- Write floats ascii or binary.
    //  If binary optionally in-place swaps argument
    void write(std::ostream& os, const bool binary, List<floatScalar>& fField);

    //- Write floats ascii or binary.
    //  If binary optionally in-place swaps argument
    void write(std::ostream&, const bool, DynamicList<floatScalar>&);

    //- Write floats ascii or binary.
    //  If binary optionally in-place swaps argument
    void write(std::ostream&, const bool, labelList&);

    //- Write floats ascii or binary.
    //  If binary optionally in-place swaps argument
    void write(std::ostream&, const bool, DynamicList<label>&);


    //- Append scalar to given DynamicList
    void insert(const scalar, DynamicList<floatScalar>&);

    //- Append point to given DynamicList
    void insert(const point&, DynamicList<floatScalar>&);

    //- Append sphericalTensor to given DynamicList
    void insert(const sphericalTensor&, DynamicList<floatScalar>&);

    //- Append symmTensor to given DynamicList
    void insert(const symmTensor&, DynamicList<floatScalar>&);

    //- Append tensor to given DynamicList
    void insert(const tensor&, DynamicList<floatScalar>&);

    //- Append elements to DynamicList
    void insert(const labelList&, DynamicList<label>&);

    //- Append elements to DynamicList
    template<class Type>
    void insert(const List<Type>&, DynamicList<floatScalar>&);

    //- Append elements of scalarList to given DynamicList using map
    void insert
    (
        const labelList& map,
        const List<scalar>& source,
        DynamicList<floatScalar>&
    );

    //- Append points to given DynamicList of floats
    void insert(const List<point>& source, DynamicList<floatScalar>&);

    //- Append points to given DynamicList of floats using map
    void insert
    (
        const labelList& map,
        const List<point>& source,
        DynamicList<floatScalar>&
    );

} // End namespace vtkWriteOps


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "vtkWriteOpsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
