/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Triangulation>
inline const Foam::Time& Foam::DelaunayMesh<Triangulation>::time() const
{
    return runTime_;
}


template<class Triangulation>
void Foam::DelaunayMesh<Triangulation>::timeCheck
(
    const string& description,
    const bool check
) const
{
    if (check)
    {
        Info<< nl << "--- [ cpuTime "
            << time().elapsedCpuTime() << " s, "
            << "delta " << time().cpuTimeIncrement()<< " s";

        if (description != word::null)
        {
            Info<< ", " << description << " ";
        }
        else
        {
            Info<< " ";
        }

        Info<< "] --- " << endl;
    }
}


template<class Triangulation>
inline Foam::label Foam::DelaunayMesh<Triangulation>::getNewCellIndex() const
{
    label id = cellCount_++;

    if (id == labelMax)
    {
        WarningInFunction
            << "Cell counter has overflowed." << endl;
    }

    return id;
}


template<class Triangulation>
inline Foam::label Foam::DelaunayMesh<Triangulation>::getNewVertexIndex() const
{
    label id = vertexCount_++;

    if (id == labelMax)
    {
        WarningInFunction
            << "Vertex counter has overflowed." << endl;
    }

    return id;
}


template<class Triangulation>
Foam::label Foam::DelaunayMesh<Triangulation>::cellCount() const
{
    return cellCount_;
}


template<class Triangulation>
Foam::label Foam::DelaunayMesh<Triangulation>::vertexCount() const
{
    return vertexCount_;
}


template<class Triangulation>
void Foam::DelaunayMesh<Triangulation>::resetCellCount()
{
    cellCount_ = 0;
}


template<class Triangulation>
void Foam::DelaunayMesh<Triangulation>::resetVertexCount()
{
    vertexCount_ = 0;
}


// ************************************************************************* //
