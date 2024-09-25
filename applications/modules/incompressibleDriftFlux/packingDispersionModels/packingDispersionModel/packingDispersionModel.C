/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

#include "none_packingDispersionModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(packingDispersionModel, 0);
    defineRunTimeSelectionTable(packingDispersionModel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::packingDispersionModel::packingDispersionModel
(
    const relativeVelocityModel& relativeVelocity
)
:
    relativeVelocity_(relativeVelocity),
    mixture_(relativeVelocity.mixture())
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::packingDispersionModel> Foam::packingDispersionModel::New
(
    const dictionary& dict,
    const relativeVelocityModel& relativeVelocity
)
{
    if (dict.found(typeName))
    {
        const word modelType(dict.lookup(typeName));

        Info<< "Selecting packing dispersion model " << modelType << endl;

        dictionaryConstructorTable::iterator cstrIter =
            dictionaryConstructorTablePtr_->find(modelType);

        if (cstrIter == dictionaryConstructorTablePtr_->end())
        {
            FatalIOErrorInFunction(dict)
                << "Unknown time scale model type " << modelType
                << ", constructor not in hash table" << nl << nl
                << "    Valid time scale model types are:" << nl
                << dictionaryConstructorTablePtr_->sortedToc()
                << abort(FatalIOError);
        }

        return
            autoPtr<packingDispersionModel>
            (
                cstrIter()
                (
                    dict.optionalSubDict(modelType + "Coeffs"),
                    relativeVelocity
                )
            );
    }
    else
    {
        return
            autoPtr<packingDispersionModel>
            (
                new packingDispersionModels::none(relativeVelocity)
            );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::packingDispersionModel::~packingDispersionModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::packingDispersionModel::Dd() const
{
    const dimensionedScalar& rhod = mixture_.rhod();
    const dimensionedScalar& rhoc = mixture_.rhoc();

    return
        rhod*mixture_.rho()
       *relativeVelocity_.UdmCoeff()*sigmaPrime()
       /((rhod - rhoc)*mixture_.alphac()*rhoc);
}


// ************************************************************************* //
