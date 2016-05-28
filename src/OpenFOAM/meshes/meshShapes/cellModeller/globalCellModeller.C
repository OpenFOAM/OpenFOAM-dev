/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

Description
    cellModeller global initializations

\*---------------------------------------------------------------------------*/

#include "cellModeller.H"
#include "etcFiles.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * Static data * * * * * * * * * * * * * * * * //

Foam::PtrList<Foam::cellModel> Foam::cellModeller::models_
(
    IFstream(findEtcFile("cellModels", true))()
);

Foam::List<Foam::cellModel*> Foam::cellModeller::modelPtrs_;

Foam::HashTable<const Foam::cellModel*> Foam::cellModeller::modelDictionary_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Construct a dummy cellModeller which reads the models and fills
    // the above tables
    cellModeller globalCellModeller_;
}

// ************************************************************************* //
