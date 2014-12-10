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

\*---------------------------------------------------------------------------*/

#include "csvSetWriter.H"
#include "coordSet.H"
#include "fileName.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::csvSetWriter<Type>::csvSetWriter()
:
    writer<Type>()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::csvSetWriter<Type>::~csvSetWriter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::fileName Foam::csvSetWriter<Type>::getFileName
(
    const coordSet& points,
    const wordList& valueSetNames
) const
{
    return this->getBaseName(points, valueSetNames) + ".csv";
}


template<class Type>
void Foam::csvSetWriter<Type>::write
(
    const coordSet& points,
    const wordList& valueSetNames,
    const List<const Field<Type>*>& valueSets,
    Ostream& os
) const
{
    writeHeader(points,valueSetNames,os);

    // Collect sets into columns
    List<const List<Type>*> columns(valueSets.size());

    forAll(valueSets, i)
    {
        columns[i] = valueSets[i];
    }

    this->writeTable(points, columns, os);
}


template<class Type>
void Foam::csvSetWriter<Type>::write
(
    const bool writeTracks,
    const PtrList<coordSet>& points,
    const wordList& valueSetNames,
    const List<List<Field<Type> > >& valueSets,
    Ostream& os
) const
{
    writeHeader(points[0],valueSetNames,os);

    if (valueSets.size() != valueSetNames.size())
    {
        FatalErrorIn("csvSetWriter<Type>::write(..)")
            << "Number of variables:" << valueSetNames.size() << endl
            << "Number of valueSets:" << valueSets.size()
            << exit(FatalError);
    }

    List<const List<Type>*> columns(valueSets.size());

    forAll(points, trackI)
    {
        // Collect sets into columns
        forAll(valueSets, i)
        {
            columns[i] = &valueSets[i][trackI];
        }

        this->writeTable(points[trackI], columns, os);
        os  << nl << nl;
    }
}


template<class Type>
void Foam::csvSetWriter<Type>::writeSeparator(Ostream& os) const
{
    os << token::COMMA;
}


namespace Foam
{
    // otherwise compiler complains about specialization
    template<>
    void csvSetWriter<scalar>::writeHeader
    (
        const coordSet& points,
        const wordList& valueSetNames,
        Ostream& os
    ) const
    {
        writeCoordHeader(points, os);

        forAll(valueSetNames, i)
        {
            if (i > 0)
            {
                writeSeparator(os);
            }
            os << valueSetNames[i];
        }

        os << nl;
    }
} // end namespace


template<class Type>
void Foam::csvSetWriter<Type>::writeHeader
(
    const coordSet& points,
    const wordList& valueSetNames,
    Ostream& os
) const
{
    writeCoordHeader(points, os);

    forAll(valueSetNames, i)
    {
        for (label j=0; j<Type::nComponents; j++)
        {
            if (i>0 || j>0)
            {
                writeSeparator(os);
            }
            os << valueSetNames[i] << "_" << j;
        }
    }

    os << nl;
}


template<class Type>
void Foam::csvSetWriter<Type>::writeCoordHeader
(
    const coordSet& points,
    Ostream& os
) const
{
    if (points.hasVectorAxis())
    {
        forAll(points, i)
        {
            os << points.axis()[i];
            writeSeparator(os);
        }
    }
    else
    {
        os << points.axis();
        writeSeparator(os);
    }
}


// ************************************************************************* //
