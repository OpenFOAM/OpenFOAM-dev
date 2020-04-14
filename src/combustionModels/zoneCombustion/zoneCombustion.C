/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2020 OpenFOAM Foundation
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

#include "zoneCombustion.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::tmp<Foam::fvScalarMatrix>
Foam::combustionModels::zoneCombustion<ReactionThermo>::filter
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


template<class ReactionThermo>
Foam::tmp<Foam::volScalarField>
Foam::combustionModels::zoneCombustion<ReactionThermo>::filter
(
    const tmp<volScalarField>& tS
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

template<class ReactionThermo>
Foam::combustionModels::zoneCombustion<ReactionThermo>::zoneCombustion
(
    const word& modelType,
    const ReactionThermo& thermo,
    const compressibleMomentumTransportModel& turb,
    const word& combustionProperties
)
:
    CombustionModel<ReactionThermo>
    (
        modelType,
        thermo,
        turb,
        combustionProperties
    ),
    combustionModelPtr_
    (
        CombustionModel<ReactionThermo>::New
        (
            thermo,
            turb,
            "zoneCombustionProperties"
        )
    ),
    zoneNames_(this->coeffs().lookup("zones"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::combustionModels::zoneCombustion<ReactionThermo>::~zoneCombustion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ReactionThermo>
const ReactionThermo&
Foam::combustionModels::zoneCombustion<ReactionThermo>::thermo() const
{
    return combustionModelPtr_->thermo();
}


template<class ReactionThermo>
void Foam::combustionModels::zoneCombustion<ReactionThermo>::correct()
{
    combustionModelPtr_->correct();
}


template<class ReactionThermo>
Foam::tmp<Foam::fvScalarMatrix>
Foam::combustionModels::zoneCombustion<ReactionThermo>::R
(
    volScalarField& Y
) const
{
    return filter(combustionModelPtr_->R(Y));
}


template<class ReactionThermo>
Foam::tmp<Foam::volScalarField>
Foam::combustionModels::zoneCombustion<ReactionThermo>::Qdot() const
{
    return filter(combustionModelPtr_->Qdot());
}


template<class ReactionThermo>
bool Foam::combustionModels::zoneCombustion<ReactionThermo>::read()
{
    if (CombustionModel<ReactionThermo>::read())
    {
        combustionModelPtr_->read();
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
