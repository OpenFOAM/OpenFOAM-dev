/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "vanDriestDelta.H"
#include "LESModel.H"
#include "wallFvPatch.H"
#include "wallDistData.H"
#include "wallPointYPlus.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(vanDriestDelta, 0);
addToRunTimeSelectionTable(LESdelta, vanDriestDelta, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void vanDriestDelta::calcDelta()
{
    const LESModel& lesModel = mesh_.lookupObject<LESModel>("LESProperties");

    const volVectorField& U = lesModel.U();
    const volScalarField& rho = lesModel.rho();
    const volScalarField& mu = lesModel.mu();
    tmp<volScalarField> muSgs = lesModel.muSgs();

    volScalarField ystar
    (
        IOobject
        (
            "ystar",
            mesh_.time().constant(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("ystar", dimLength, GREAT)
    );

    const fvPatchList& patches = mesh_.boundary();
    forAll(patches, patchi)
    {
        if (isA<wallFvPatch>(patches[patchi]))
        {
            const fvPatchVectorField& Uw = U.boundaryField()[patchi];
            const scalarField& rhow = rho.boundaryField()[patchi];
            const scalarField& muw = mu.boundaryField()[patchi];
            const scalarField& muSgsw = muSgs().boundaryField()[patchi];

            ystar.boundaryField()[patchi] =
                muw/(rhow*sqrt((muw + muSgsw)*mag(Uw.snGrad())/rhow + VSMALL));
        }
    }

    scalar cutOff = wallPointYPlus::yPlusCutOff;
    wallPointYPlus::yPlusCutOff = 500;
    wallDistData<wallPointYPlus> y(mesh_, ystar);
    wallPointYPlus::yPlusCutOff = cutOff;

    delta_ = min
    (
        static_cast<const volScalarField&>(geometricDelta_()),
        (kappa_/Cdelta_)*((scalar(1) + SMALL) - exp(-y/ystar/Aplus_))*y
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

vanDriestDelta::vanDriestDelta
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dd
)
:
    LESdelta(name, mesh),
    geometricDelta_
    (
        LESdelta::New("geometricDelta", mesh, dd.subDict(type() + "Coeffs"))
    ),
    kappa_(dd.lookupOrDefault<scalar>("kappa", 0.41)),
    Aplus_
    (
        dd.subDict(type() + "Coeffs").lookupOrDefault<scalar>("Aplus", 26.0)
    ),
    Cdelta_
    (
        dd.subDict(type() + "Coeffs").lookupOrDefault<scalar>("Cdelta", 0.158)
    ),
    calcInterval_
    (
        dd.subDict(type() + "Coeffs").lookupOrDefault<label>("calcInterval", 1)
    )
{
    delta_ = geometricDelta_();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void vanDriestDelta::read(const dictionary& d)
{
    const dictionary& dd(d.subDict(type() + "Coeffs"));

    geometricDelta_().read(dd);
    d.readIfPresent<scalar>("kappa", kappa_);
    dd.readIfPresent<scalar>("Aplus", Aplus_);
    dd.readIfPresent<scalar>("Cdelta", Cdelta_);
    dd.readIfPresent<label>("calcInterval", calcInterval_);
    calcDelta();
}


void vanDriestDelta::correct()
{
    if (mesh().time().timeIndex() % calcInterval_ == 0)
    {
        geometricDelta_().correct();
        calcDelta();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
