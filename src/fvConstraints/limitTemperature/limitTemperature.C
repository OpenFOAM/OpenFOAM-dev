/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2021 OpenFOAM Foundation
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
        fvConstraint,
        limitTemperature,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::limitTemperature::readCoeffs()
{
    Tmin_ = coeffs().lookup<scalar>("min");
    Tmax_ = coeffs().lookup<scalar>("max");
    phaseName_ = coeffs().lookupOrDefault<word>("phase", word::null);
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
    fvConstraint(name, modelType, dict, mesh),
    set_(coeffs(), mesh),
    Tmin_(-vGreat),
    Tmax_(vGreat),
    phaseName_(word::null)
{
    readCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::limitTemperature::constrainedFields() const
{
    const basicThermo& thermo =
        mesh().lookupObject<basicThermo>
        (
            IOobject::groupName(basicThermo::dictName, phaseName_)
        );

    return wordList(1, thermo.he().name());
}


bool Foam::fv::limitTemperature::constrain(volScalarField& he) const
{
    const basicThermo& thermo =
        mesh().lookupObject<basicThermo>
        (
            IOobject::groupName(basicThermo::dictName, phaseName_)
        );

    const labelList& cells = set_.cells();

    scalarField Tmin(cells.size(), Tmin_);
    scalarField Tmax(cells.size(), Tmax_);

    scalarField heMin(thermo.he(Tmin, cells));
    scalarField heMax(thermo.he(Tmax, cells));

    scalarField& hec = he.primitiveFieldRef();

    forAll(cells, i)
    {
        const label celli = cells[i];
        hec[celli]= max(min(hec[celli], heMax[i]), heMin[i]);
    }

    // Handle boundaries in the case of 'all'
    if (set_.selectionMode() == fvCellSet::selectionModeType::all)
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

    return cells.size();
}


void Foam::fv::limitTemperature::updateMesh(const mapPolyMesh& mpm)
{
    set_.updateMesh(mpm);
}


bool Foam::fv::limitTemperature::read(const dictionary& dict)
{
    if (fvConstraint::read(dict))
    {
        set_.read(coeffs());
        readCoeffs();
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
