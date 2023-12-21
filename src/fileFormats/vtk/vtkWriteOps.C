/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "vtkWriteOps.H"

#if defined(__mips)
    #include <standards.h>
    #include <sys/endian.h>
#endif

// MacOSX
#ifdef __DARWIN_BYTE_ORDER
    #if __DARWIN_BYTE_ORDER==__DARWIN_BIG_ENDIAN
        #undef LITTLE_ENDIAN
    #else
        #undef BIG_ENDIAN
    #endif
#endif

#if defined(LITTLE_ENDIAN) \
 || defined(_LITTLE_ENDIAN) \
 || defined(__LITTLE_ENDIAN)
    #define LITTLEENDIAN 1
#elif defined(BIG_ENDIAN) || defined(_BIG_ENDIAN) || defined(__BIG_ENDIAN)
    #undef LITTLEENDIAN
#else
    #error "Cannot find LITTLE_ENDIAN or BIG_ENDIAN symbol defined."
    #error "Please add to compilation options"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::vtkWriteOps::swapWord(label& word32)
{
    char* mem = reinterpret_cast<char*>(&word32);

    char a = mem[0];
    mem[0] = mem[3];
    mem[3] = a;

    a = mem[1];
    mem[1] = mem[2];
    mem[2] = a;
}


void Foam::vtkWriteOps::swapWords(const label nWords, label* words32)
{
    for (label i = 0; i < nWords; i++)
    {
        swapWord(words32[i]);
    }
}


void Foam::vtkWriteOps::write
(
    std::ostream& os,
    const bool binary,
    List<floatScalar>& fField
)
{
    if (binary)
    {
        #ifdef LITTLEENDIAN
        swapWords(fField.size(), reinterpret_cast<label*>(fField.begin()));
        #endif
        os.write
        (
            reinterpret_cast<char*>(fField.begin()),
            fField.size()*sizeof(floatScalar)
        );

        os  << std::endl;
    }
    else
    {
        forAll(fField, i)
        {
            os  << fField[i];

            if (i > 0 && (i % 10) == 0)
            {
                os  << std::endl;
            }
            else
            {
                os  << ' ';
            }
        }
        os << std::endl;
    }
}


void Foam::vtkWriteOps::write
(
    std::ostream& os,
    const bool binary,
    DynamicList<floatScalar>& fField
)
{
    List<floatScalar>& fld = fField.shrink();

    write(os, binary, fld);
}


void Foam::vtkWriteOps::write
(
    std::ostream& os,
    const bool binary,
    labelList& elems
)
{
    if (binary)
    {
        #ifdef LITTLEENDIAN
        swapWords(elems.size(), reinterpret_cast<label*>(elems.begin()));
        #endif
        os.write
        (
            reinterpret_cast<char*>(elems.begin()),
            elems.size()*sizeof(label)
        );

        os  << std::endl;
    }
    else
    {
        forAll(elems, i)
        {
            os  << elems[i];

            if (i > 0 && (i % 10) == 0)
            {
                os  << std::endl;
            }
            else
            {
                os  << ' ';
            }
        }
        os  << std::endl;
    }
}


void Foam::vtkWriteOps::write
(
    std::ostream& os,
    const bool binary,
    DynamicList<label>& elems
)
{
    labelList& fld = elems.shrink();

    write(os, binary, fld);
}


void Foam::vtkWriteOps::writeHeader
(
    std::ostream& os,
    const bool binary,
    const std::string& title
)
{
    os  << "# vtk DataFile Version 2.0" << std::endl
        << title << std::endl;

    if (binary)
    {
        os  << "BINARY" << std::endl;
    }
    else
    {
        os  << "ASCII" << std::endl;
    }
}


void Foam::vtkWriteOps::writeCellDataHeader
(
    std::ostream& os,
    const label nCells,
    const label nFields
)
{
    os  << "CELL_DATA " << nCells << std::endl
        << "FIELD attributes " << nFields << std::endl;
}


void Foam::vtkWriteOps::writePointDataHeader
(
    std::ostream& os,
    const label nPoints,
    const label nFields
)
{
    os  << "POINT_DATA  " << nPoints << std::endl
        << "FIELD attributes " << nFields << std::endl;
}


void Foam::vtkWriteOps::insert(const scalar src, DynamicList<floatScalar>& dest)
{
    dest.append(cast(src));
}


void Foam::vtkWriteOps::insert
(
    const vector& src,
    DynamicList<floatScalar>& dest
)
{
    for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
    {
        dest.append(cast(src[cmpt]));
    }
}


void Foam::vtkWriteOps::insert
(
    const sphericalTensor& src,
    DynamicList<floatScalar>& dest
)
{
    for (direction cmpt = 0; cmpt < sphericalTensor::nComponents; ++cmpt)
    {
        dest.append(cast(src[cmpt]));
    }
}


void Foam::vtkWriteOps::insert
(
    const symmTensor& src,
    DynamicList<floatScalar>& dest
)
{
    dest.append(cast(src.xx()));
    dest.append(cast(src.yy()));
    dest.append(cast(src.zz()));
    dest.append(cast(src.xy()));
    dest.append(cast(src.yz()));
    dest.append(cast(src.xz()));
}


void Foam::vtkWriteOps::insert
(
    const tensor& src,
    DynamicList<floatScalar>& dest
)
{
    for (direction cmpt = 0; cmpt < tensor::nComponents; ++cmpt)
    {
        dest.append(cast(src[cmpt]));
    }
}


void Foam::vtkWriteOps::insert(const labelList& src, DynamicList<label>& dest)
{
    dest.append(src);
}


void Foam::vtkWriteOps::insert
(
    const labelList& map,
    const List<scalar>& source,
    DynamicList<floatScalar>& dest
)
{
    forAll(map, i)
    {
        dest.append(cast(source[map[i]]));
    }
}


void Foam::vtkWriteOps::insert
(
    const List<point>& source,
    DynamicList<floatScalar>& dest
)
{
    forAll(source, i)
    {
       insert(source[i], dest);
    }
}


void Foam::vtkWriteOps::insert
(
    const labelList& map,
    const List<point>& source,
    DynamicList<floatScalar>& dest
)
{
    forAll(map, i)
    {
       insert(source[map[i]], dest);
    }
}


// ************************************************************************* //
