/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "CloudFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::autoPtr<Foam::CloudFunctionObject<CloudType> >
Foam::CloudFunctionObject<CloudType>::New
(
    const dictionary& dict,
    CloudType& owner,
    const word& objectType,
    const word& modelName
)
{
    Info<< "    Selecting cloud function " << modelName << " of type "
        << objectType << endl;

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(objectType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "CloudFunctionObject<CloudType>::New"
            "("
                "const dictionary&, "
                "CloudType&, "
                "const word&, "
                "const word&"
            ")"
        )   << "Unknown cloud function type "
            << objectType << nl << nl
            << "Valid cloud function types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<CloudFunctionObject<CloudType> >
    (
        cstrIter()
        (
            dict,
            owner,
            modelName
        )
    );
}


// ************************************************************************* //
