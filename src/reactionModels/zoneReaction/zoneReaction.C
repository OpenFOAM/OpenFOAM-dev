/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2026 OpenFOAM Foundation
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

#include "zoneReaction.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace reactionModels
{
    defineTypeNameAndDebug(zoneReaction, 0);
    addToRunTimeSelectionTable(reactionModel, zoneReaction, dictionary);
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::fvScalarMatrix>
Foam::reactionModels::zoneReaction::filter
(
    const tmp<fvScalarMatrix>& tR
) const
{
    fvScalarMatrix& R = tR.ref();
    scalarField& Su = R.source();
    scalarField filteredField(Su.size(), 0);

    forAll(zoneNames_, zonei)
    {
        const labelList& cells = this->mesh().cellZones()[zoneNames_[zonei]];

        forAll(cells, i)
        {
            filteredField[cells[i]] = Su[cells[i]];
        }
    }

    Su = filteredField;

    if (R.hasDiag())
    {
        scalarField& Sp = R.diag();

        forAll(zoneNames_, zonei)
        {
            const labelList& cells =
                this->mesh().cellZones()[zoneNames_[zonei]];

            forAll(cells, i)
            {
                filteredField[cells[i]] = Sp[cells[i]];
            }
        }

        Sp = filteredField;
    }

    return tR;
}


template<class GeoField>
inline Foam::tmp<GeoField>
Foam::reactionModels::zoneReaction::filter
(
    const tmp<GeoField>& tS
) const
{
    scalarField& S = tS.ref();
    scalarField filteredField(S.size(), 0);

    forAll(zoneNames_, zonei)
    {
        const labelList& cells = this->mesh().cellZones()[zoneNames_[zonei]];

        forAll(cells, i)
        {
            filteredField[cells[i]] = S[cells[i]];
        }
    }

    S = filteredField;

    return tS;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reactionModels::zoneReaction::zoneReaction
(
    const word& modelType,
    const fluidMulticomponentThermo& thermo,
    const compressibleMomentumTransportModel& turb,
    const word& reactionProperties
)
:
    reactionModel
    (
        modelType,
        thermo,
        turb,
        reactionProperties
    ),
    reactionModelPtr_
    (
        reactionModel::New
        (
            thermo,
            turb,
            "zoneReactionProperties"
        )
    ),
    zoneNames_(this->coeffs().lookup("zones"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::reactionModels::zoneReaction::~zoneReaction()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::reactionModels::zoneReaction::correct()
{
    reactionModelPtr_->correct();
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::reactionModels::zoneReaction::R(const label speciei) const
{
    return filter(reactionModelPtr_->R(speciei));
}


Foam::tmp<Foam::fvScalarMatrix>
Foam::reactionModels::zoneReaction::R(volScalarField& Y) const
{
    return filter(reactionModelPtr_->R(Y));
}


Foam::tmp<Foam::volScalarField>
Foam::reactionModels::zoneReaction::Qdot() const
{
    return filter(reactionModelPtr_->Qdot());
}


bool Foam::reactionModels::zoneReaction::read()
{
    if (reactionModel::read())
    {
        reactionModelPtr_->read();
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
