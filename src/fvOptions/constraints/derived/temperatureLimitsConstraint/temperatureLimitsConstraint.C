/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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

#include "temperatureLimitsConstraint.H"
#include "fvMesh.H"
#include "basicThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(temperatureLimitsConstraint, 0);
    addToRunTimeSelectionTable
    (
        option,
        temperatureLimitsConstraint,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::temperatureLimitsConstraint::temperatureLimitsConstraint
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    option(name, modelType, dict, mesh),
    Tmin_(readScalar(coeffs_.lookup("Tmin"))),
    Tmax_(readScalar(coeffs_.lookup("Tmax")))
{
    fieldNames_.setSize(1, "energy");
    applied_.setSize(1, false);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::temperatureLimitsConstraint::alwaysApply() const
{
    return true;
}


void Foam::fv::temperatureLimitsConstraint::correct(volScalarField& he)
{
    const basicThermo& thermo =
        mesh_.lookupObject<basicThermo>("thermophysicalProperties");

    if (he.name() == thermo.he().name())
    {
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

            forAll(bf, patchI)
            {
                fvPatchScalarField& hep = bf[patchI];
                if (hep.fixesValue())
                {
                    // not over-riding fixed conditions
                    continue;
                }

                const scalarField& pp = thermo.p().boundaryField()[patchI];

                scalarField Tminp(pp.size(), Tmin_);
                scalarField Tmaxp(pp.size(), Tmax_);

                scalarField heMinp(thermo.he(pp, Tminp, patchI));
                scalarField heMaxp(thermo.he(pp, Tmaxp, patchI));

                forAll(hep, faceI)
                {
                    hep[faceI] =
                        max(min(hep[faceI], heMaxp[faceI]), heMinp[faceI]);
                }
            }
        }
    }
}


void Foam::fv::temperatureLimitsConstraint::writeData(Ostream& os) const
{
    os  << indent << name_ << endl;
    dict_.write(os);
}


bool Foam::fv::temperatureLimitsConstraint::read(const dictionary& dict)
{
    if (option::read(dict))
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
