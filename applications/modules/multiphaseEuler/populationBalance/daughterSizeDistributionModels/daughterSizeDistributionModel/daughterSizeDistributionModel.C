/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2025 OpenFOAM Foundation
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

#include "daughterSizeDistributionModel.H"
#include "breakupModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalance
{
    defineTypeNameAndDebug(daughterSizeDistributionModel, 0);
    defineRunTimeSelectionTable(daughterSizeDistributionModel, dictionary);
}
}


// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::populationBalance::daughterSizeDistributionModel>
Foam::populationBalance::daughterSizeDistributionModel::New
(
    const breakupModel& breakup,
    const dictionary& dict
)
{
    word daughterSizeDistributionModelType
    (
        dict.lookup("daughterSizeDistributionModel")
    );

    Info<< "Selecting daughter size distribution model for "
        << breakup.popBal().name() << ": "
        << daughterSizeDistributionModelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(daughterSizeDistributionModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown daughter size distribution model type "
            << daughterSizeDistributionModelType << endl << endl
            << "Valid daughter size distribution model types are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(breakup, dict);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::populationBalance::daughterSizeDistributionModel::
daughterSizeDistributionModel
(
    const breakupModel& breakup,
    const dictionary& dict
)
:
    breakup_(breakup),
    nik_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalance::daughterSizeDistributionModel::
~daughterSizeDistributionModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::dimensionedScalar&
Foam::populationBalance::daughterSizeDistributionModel::nik
(
    const label i,
    const label k
) const
{
    return nik_[k][i];
}


void Foam::populationBalance::daughterSizeDistributionModel::precompute()
{
    if (nik_.size() == 0)
    {
        forAll(breakup_.popBal().phases(), k)
        {
            nik_.append(new PtrList<dimensionedScalar>());

            for (label i = 0; i <= k; i++)
            {
                nik_[k].append
                (
                    new dimensionedScalar
                    (
                        this->calcNik(i, k)
                    )
                );
            }
        }
    }
}


// ************************************************************************* //
