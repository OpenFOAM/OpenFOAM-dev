/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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

#include "kappaEff.H"
#include "kinematicMomentumTransportModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace wallHeatTransferCoeffModels
{
    defineTypeNameAndDebug(kappaEff, 0);
    addToRunTimeSelectionTable
    (
        wallHeatTransferCoeffModel,
        kappaEff,
        word
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallHeatTransferCoeffModels::kappaEff::kappaEff
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    wallHeatTransferCoeffModel(name, mesh, dict),
    mesh_(mesh),
    Prl_("Prl", dimless, dict),
    Prt_("Prl", dimless, dict),
    Lchar_("Lchar", dimLength, Zero),
    isCharLength_(false)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wallHeatTransferCoeffModels::kappaEff::~kappaEff()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::wallHeatTransferCoeffModels::kappaEff::read(const dictionary& dict)
{
    Prl_.read(dict);
    Prt_.read(dict);
    isCharLength_ = dict.found("Lchar");

    if (isCharLength_)
    {
        Lchar_.read(dict);
    }

    return true;
}


Foam::tmp<Foam::volScalarField>
Foam::wallHeatTransferCoeffModels::kappaEff::htcByRhoCp
(
    const momentumTransportModel& mmtm,
    const labelHashSet& patches
) const
{
    tmp<volScalarField> thtcByRhoCp
    (
        volScalarField::New
        (
            type(),
            mesh_,
            isCharLength_
          ? dimensionedScalar(dimLength/dimTime, 0)
          : dimensionedScalar(dimArea/dimTime, 0)
        )
    );

    volScalarField::Boundary& thtcByRhoCpBf =
        thtcByRhoCp.ref().boundaryFieldRef();

    forAllConstIter(labelHashSet, patches, iter)
    {
        label patchi = iter.key();

        if (!thtcByRhoCpBf[patchi].coupled())
        {
            thtcByRhoCpBf[patchi] =
                mmtm.nu(patchi)/Prl_.value() + mmtm.nut(patchi)/Prt_.value();

            if (isCharLength_)
            {
                thtcByRhoCpBf[patchi] /= Lchar_.value();
            }
        }
    }
    return thtcByRhoCp;

}

// ************************************************************************* //
