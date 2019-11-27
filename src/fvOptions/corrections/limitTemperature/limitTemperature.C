/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2019 OpenFOAM Foundation
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
    Tmin_(coeffs_.lookup<scalar>("min")),
    Tmax_(coeffs_.lookup<scalar>("max")),
    phase_(coeffs_.lookupOrDefault<word>("phase", word::null))
{
    // Set the field name to that of the energy field from which the temperature
    // is obtained
    const basicThermo& thermo =
        mesh_.lookupObject<basicThermo>
        (
            IOobject::groupName(basicThermo::dictName, phase_)
        );

    fieldNames_.setSize(1, thermo.he().name());

    applied_.setSize(1, false);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::limitTemperature::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        coeffs_.lookup("min") >> Tmin_;
        coeffs_.lookup("max") >> Tmax_;

        return true;
    }
    else
    {
        return false;
    }
}


void Foam::fv::limitTemperature::correct(volScalarField& he)
{
    const basicThermo& thermo =
        mesh_.lookupObject<basicThermo>
        (
            IOobject::groupName(basicThermo::dictName, phase_)
        );

    scalarField Tmin(cells_.size(), Tmin_);
    scalarField Tmax(cells_.size(), Tmax_);

    scalarField heMin(thermo.he(Tmin, cells_));
    scalarField heMax(thermo.he(Tmax, cells_));

    scalarField& hec = he.primitiveFieldRef();

    forAll(cells_, i)
    {
        label celli = cells_[i];
        hec[celli]= max(min(hec[celli], heMax[i]), heMin[i]);
    }

    // handle boundaries in the case of 'all'
    if (selectionMode_ == smAll)
    {
        volScalarField::Boundary& bf = he.boundaryFieldRef();

        forAll(bf, patchi)
        {
            fvPatchScalarField& hep = bf[patchi];

            if (!hep.fixesValue())
            {
                scalarField Tminp(hep.size(), Tmin_);
                scalarField Tmaxp(hep.size(), Tmax_);

                scalarField heMinp(thermo.he(Tminp, patchi));
                scalarField heMaxp(thermo.he(Tmaxp, patchi));

                forAll(hep, facei)
                {
                    hep[facei] =
                        max(min(hep[facei], heMaxp[facei]), heMinp[facei]);
                }
            }
        }
    }
}


// ************************************************************************* //
