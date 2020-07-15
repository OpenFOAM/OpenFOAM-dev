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

#include "IATEphaseChange.H"
#include "phaseSystem.H"
#include "ThermalPhaseChangePhaseSystem.H"
#include "MomentumTransferPhaseSystem.H"
#include "fvmSup.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace IATEsources
{
    defineTypeNameAndDebug(phaseChange, 0);
    addToRunTimeSelectionTable(IATEsource, phaseChange, dictionary);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::IATEsources::phaseChange::phaseChange
(
    const IATE& iate,
    const dictionary& dict
)
:
    IATEsource(iate),
    otherPhaseName_(dict.lookup("otherPhase")),
    dmdtfName_(dict.lookup("dmdtf")),
    specieName_(dict.lookupOrDefault<word>("specie", word::null))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::fvScalarMatrix>
Foam::diameterModels::IATEsources::phaseChange::R
(
    const volScalarField& alphai,
    volScalarField& kappai
) const
{
    const phasePair& pair =
        phase().fluid().phasePairs()
        [
            phasePairKey(phase().name(), otherPhaseName_)
        ];

    const volScalarField& dmdtf =
        alphai.mesh().lookupObject<volScalarField>
        (
            IOobject::groupName
            (
                IOobject::groupName(dmdtfName_, specieName_),
                pair.name()
            )
        );

    const label dmdtfSign =
        phase().name() == pair.first() ? +1 : -1;

    return -fvm::SuSp(dmdtfSign*dmdtf/(3*alphai*phase().rho()), kappai);
}


// ************************************************************************* //
