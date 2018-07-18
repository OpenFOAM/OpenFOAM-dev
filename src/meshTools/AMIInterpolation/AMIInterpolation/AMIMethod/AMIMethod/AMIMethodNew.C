/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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

#include "AMIMethod.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::AMIMethod> Foam::AMIMethod::New
(
    const word& methodName,
    const primitivePatch& srcPatch,
    const primitivePatch& tgtPatch,
    const scalarField& srcMagSf,
    const scalarField& tgtMagSf,
    const faceAreaIntersect::triangulationMode& triMode,
    const bool reverseTarget,
    const bool requireMatch
)
{
    if (debug)
    {
        Info<< "Selecting AMIMethod " << methodName << endl;
    }

    typename componentsConstructorTable::iterator cstrIter =
        componentsConstructorTablePtr_->find(methodName);

    if (cstrIter == componentsConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown AMIMethod type "
            << methodName << nl << nl
            << "Valid AMIMethod types are:" << nl
            << componentsConstructorTablePtr_->sortedToc() << exit(FatalError);
    }

    return autoPtr<AMIMethod>
    (
        cstrIter()
        (
            srcPatch,
            tgtPatch,
            srcMagSf,
            tgtMagSf,
            triMode,
            reverseTarget,
            requireMatch
        )
    );
}


// ************************************************************************* //
