/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2015 OpenFOAM Foundation
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

#include "limitTemperature.H"
#include "fvMesh.H"
#include "basicThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(limitTemperature, 0);
    addToRunTimeSelectionTable
    (
        option,
        limitTemperature,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::limitTemperature::limitTemperature
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    cellSetOption(name, modelType, dict, mesh),
    Tmin_(readScalar(coeffs_.lookup("Tmin"))),
    Tmax_(readScalar(coeffs_.lookup("Tmax")))
{
    // Set the field name to that of the energy field from which the temperature
    // is obtained

    const basicThermo& thermo =
        mesh_.lookupObject<basicThermo>(basicThermo::dictName);

    fieldNames_.setSize(1, thermo.he().name());

    applied_.setSize(1, false);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::limitTemperature::correct(volScalarField& he)
{
    const basicThermo& thermo =
        mesh_.lookupObject<basicThermo>(basicThermo::dictName);

    scalarField Tmin(cells_.size(), Tmin_);
    scalarField Tmax(cells_.size(), Tmax_);

    scalarField heMin(thermo.he(thermo.p(), Tmin, cells_));
    scalarField heMax(thermo.he(thermo.p(), Tmax, cells_));

    scalarField& hec = he.internalField();

    forAll(cells_, i)
    {
        label cellI = cells_[i];
        hec[cellI]= max(min(hec[cellI], heMax[i]), heMin[i]);
    }

    // handle boundaries in the case of 'all'
    if (selectionMode_ == smAll)
    {
        volScalarField::GeometricBoundaryField& bf = he.boundaryField();

        forAll(bf, patchi)
        {
            fvPatchScalarField& hep = bf[patchi];

            if (!hep.fixesValue())
            {
                const scalarField& pp = thermo.p().boundaryField()[patchi];

                scalarField Tminp(pp.size(), Tmin_);
                scalarField Tmaxp(pp.size(), Tmax_);

                scalarField heMinp(thermo.he(pp, Tminp, patchi));
                scalarField heMaxp(thermo.he(pp, Tmaxp, patchi));

                forAll(hep, facei)
                {
                    hep[facei] =
                        max(min(hep[facei], heMaxp[facei]), heMinp[facei]);
                }
            }
        }
    }
}


bool Foam::fv::limitTemperature::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        coeffs_.readIfPresent("Tmin", Tmin_);
        coeffs_.readIfPresent("Tmax", Tmax_);

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
