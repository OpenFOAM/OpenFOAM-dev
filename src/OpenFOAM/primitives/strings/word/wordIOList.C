/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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

#include "wordIOList.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineCompoundTypeName(List<word>, wordList);
    addCompoundToRunTimeSelectionTable(List<word>, wordList);

    defineTemplateTypeNameAndDebugWithName(wordIOList, "wordList", 0);
    defineTemplateTypeNameAndDebugWithName(wordListIOList, "wordListList", 0);
}


void Foam::printTable
(
    const List<wordList>& wll,
    List<string::size_type>& columnWidth,
    Ostream& os
)
{
    if (!wll.size()) return;

    // Find the maximum word length for each column
    columnWidth.setSize(wll[0].size(), string::size_type(0));
    forAll(columnWidth, j)
    {
        forAll(wll, i)
        {
            columnWidth[j] = max(columnWidth[j], wll[i][j].size());
        }
    }

    // Print the rows adding spacing for the columns
    forAll(wll, i)
    {
        forAll(wll[i], j)
        {
            os  << wll[i][j];
            for
            (
                string::size_type k=0;
                k<columnWidth[j] - wll[i][j].size() + 2;
                k++
            )
            {
                os  << ' ';
            }
        }
        os  << nl;

        if (i == 0) os  << nl;
    }
}


void Foam::printTable(const List<wordList>& wll, Ostream& os)
{
    List<string::size_type> columnWidth;
    printTable(wll, columnWidth, os);
}


// ************************************************************************* //
