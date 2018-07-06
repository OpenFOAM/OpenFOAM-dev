/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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

\*---------------------------------------------------------------------------*/

#include "gnuplotGraph.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(gnuplotGraph, 0);
    const word gnuplotGraph::ext_("gplt");

    typedef graph::writer graphWriter;
    addToRunTimeSelectionTable(graphWriter, gnuplotGraph, word);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::gnuplotGraph::write(const graph& g, Ostream& os) const
{
    os  << "#set term postscript color" << endl
        << "set output \"" << word(g.title()) << ".ps\"" << endl
        << "set title " << g.title() << " 0,0" << endl << "show title" << endl
        << "set xlabel " << g.xName() << " 0,0" << endl << "show xlabel" << endl
        << "set ylabel " << g.yName() << " 0,0" << endl << "show ylabel" << endl
        << "plot";

    bool firstField = true;

    forAllConstIter(graph, g, iter)
    {
        if (!firstField)
        {
            os << ',';
        }
        firstField = false;

        os  << "'-' title " << iter()->name() << " with lines";
    }
    os << "; pause -1" << endl;


    forAllConstIter(graph, g, iter)
    {
        os  << endl;
        writeXY(g.x(), *iter(), os);
    }
}


// ************************************************************************* //
