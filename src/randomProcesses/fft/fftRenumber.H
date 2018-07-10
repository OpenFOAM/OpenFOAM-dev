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

Global
    fftRenumber

Description
    Multi-dimensional renumbering used in the Numerical Recipes
    fft routine.

SourceFiles
    fftRenumber.C

\*---------------------------------------------------------------------------*/

#ifndef fftRenumber_H
#define fftRenumber_H

#include "complex.H"
#include "List.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Recursively evaluate the indexing necessary to do the folding of the fft
// data. We recurse until we have the indexing ncessary for the folding in all
// directions.
void fftRenumberRecurse
(
    List<complex>& data,
    List<complex>& renumData,
    const labelList& nn,
    label nnprod,
    label ii,
    label l1,
    label l2
);


// Fold the n-d data array to get the fft components in the right places.
void fftRenumber
(
    List<complex>& data,
    const labelList& nn
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
