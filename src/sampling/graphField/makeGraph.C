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

Global
    makeGraph

Description
    Write a graph file for a field given the data point locations field,
    the field of interest and the name of the field to be used for the
    graph file name.

\*---------------------------------------------------------------------------*/

#include "makeGraph.H"
#include "volFields.H"
#include "graph.H"
#include "writeFile.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void makeGraph
(
    const scalarField& x,
    const volScalarField& vsf,
    const word& graphFormat
)
{
    makeGraph(x, vsf, vsf.name(), graphFormat);
}


void makeGraph
(
    const scalarField& x,
    const volScalarField& vsf,
    const word& name,
    const word& graphFormat
)
{
    fileName path
    (
        vsf.rootPath()
       /vsf.caseName()
       /functionObjects::writeFile::outputPrefix
       /"graphs"
       /vsf.instance()
    );
    mkDir(path);

    makeGraph
    (
        x,
        vsf.primitiveField(),
        name,
        path,
        graphFormat
    );
}


void makeGraph
(
    const scalarField& x,
    const scalarField& sf,
    const word& name,
    const fileName& path,
    const word& graphFormat
)
{
    graph
    (
        name,
        "x",
        name,
        x,
        sf
    ).write(path/name, graphFormat);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
