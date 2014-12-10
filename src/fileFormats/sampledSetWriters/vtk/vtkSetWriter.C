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

#include "vtkSetWriter.H"
#include "coordSet.H"
#include "fileName.H"
#include "OFstream.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::vtkSetWriter<Type>::vtkSetWriter()
:
    writer<Type>()
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::vtkSetWriter<Type>::~vtkSetWriter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::fileName Foam::vtkSetWriter<Type>::getFileName
(
    const coordSet& points,
    const wordList& valueSetNames
) const
{
    return this->getBaseName(points, valueSetNames) + ".vtk";
}


template<class Type>
void Foam::vtkSetWriter<Type>::write
(
    const coordSet& points,
    const wordList& valueSetNames,
    const List<const Field<Type>*>& valueSets,
    Ostream& os
) const
{
    os  << "# vtk DataFile Version 2.0" << nl
        << points.name() << nl
        << "ASCII" << nl
        << "DATASET POLYDATA" << nl
        << "POINTS " << points.size() << " float" << nl;

    forAll(points, i)
    {
        const vector& pt = points[i];
        os  << pt.x() << ' ' << pt.y() << ' ' << pt.z() << nl;
    }


    os  << "POINT_DATA " << points.size() << nl
        << " FIELD attributes " << valueSetNames.size() << nl;

    forAll(valueSetNames, setI)
    {
        os  << valueSetNames[setI] << ' ' << pTraits<Type>::nComponents << ' '
            << points.size() << " float" << nl;

        const Field<Type>& fld = *valueSets[setI];

        forAll(fld, pointI)
        {
            if (pointI != 0)
            {
                os  << ' ';
            }
            writer<Type>::write(fld[pointI], os);
        }
        os  << nl;
    }
}


template<class Type>
void Foam::vtkSetWriter<Type>::write
(
    const bool writeTracks,
    const PtrList<coordSet>& tracks,
    const wordList& valueSetNames,
    const List<List<Field<Type> > >& valueSets,
    Ostream& os
) const
{
    if (valueSets.size() != valueSetNames.size())
    {
        FatalErrorIn("vtkSetWriter<Type>::write(..)")
            << "Number of variables:" << valueSetNames.size() << endl
            << "Number of valueSets:" << valueSets.size()
            << exit(FatalError);
    }

    label nTracks = tracks.size();
    label nPoints = 0;
    forAll(tracks, i)
    {
        nPoints += tracks[i].size();
    }

    os  << "# vtk DataFile Version 2.0" << nl
        << tracks[0].name() << nl
        << "ASCII" << nl
        << "DATASET POLYDATA" << nl
        << "POINTS " << nPoints << " float" << nl;

    forAll(tracks, trackI)
    {
        const coordSet& points = tracks[trackI];
        forAll(points, i)
        {
            const vector& pt = points[i];
            os  << pt.x() << ' ' << pt.y() << ' ' << pt.z() << nl;
        }
    }

    if (writeTracks)
    {
        os  << "LINES " << nTracks << ' ' << nPoints+nTracks << nl;

        // Write ids of track points to file
        label globalPtI = 0;
        forAll(tracks, trackI)
        {
            const coordSet& points = tracks[trackI];

            os  << points.size();
            forAll(points, i)
            {
                os  << ' ' << globalPtI;
                globalPtI++;
            }
            os << nl;
        }
    }

    os  << "POINT_DATA " << nPoints << nl
        << " FIELD attributes " << valueSetNames.size() << nl;

    forAll(valueSetNames, setI)
    {
        os  << valueSetNames[setI] << ' ' << pTraits<Type>::nComponents << ' '
            << nPoints << " float" << nl;

        const List<Field<Type> >& fieldVals = valueSets[setI];

        forAll(fieldVals, i)
        {
            const Field<Type>& vals = fieldVals[i];

            forAll(vals, j)
            {
                if (j != 0)
                {
                    os  << ' ';
                }
                writer<Type>::write(vals[j], os);
            }
            os  << nl;
        }
    }
}


// ************************************************************************* //
