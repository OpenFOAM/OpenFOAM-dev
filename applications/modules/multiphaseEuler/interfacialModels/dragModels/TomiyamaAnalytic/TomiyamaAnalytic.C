/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2023 OpenFOAM Foundation
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

#include "TomiyamaAnalytic.H"
#include "aspectRatioModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(TomiyamaAnalytic, 0);
    addToRunTimeSelectionTable(dragModel, TomiyamaAnalytic, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::TomiyamaAnalytic::TomiyamaAnalytic
(
    const dictionary& dict,
    const phaseInterface& interface,
    const bool registerObject
)
:
    dispersedDragModel(dict, interface, registerObject),
    residualEo_("residualEo", dimless, dict),
    residualE_("residualE", dimless, dict),
    aspectRatio_(aspectRatioModel::New(dict.subDict("aspectRatio"), interface))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dragModels::TomiyamaAnalytic::~TomiyamaAnalytic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::dragModels::TomiyamaAnalytic::CdRe() const
{
    const volScalarField Eo(max(interface_.Eo(), residualEo_));
    const volScalarField E(max(aspectRatio_->E(), residualE_));

    const volScalarField OmEsq(max(1 - sqr(E), sqr(residualE_)));
    const volScalarField rtOmEsq(sqrt(OmEsq));

    const volScalarField F(max(asin(rtOmEsq) - E*rtOmEsq, residualE_)/OmEsq);

    return
        (8.0/3.0)
       *Eo/(Eo*pow(E, 2.0/3.0)/OmEsq + 16*pow(E, 4.0/3.0))/sqr(F)
       *interface_.Re();
}


// ************************************************************************* //
