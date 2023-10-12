/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2023 OpenFOAM Foundation
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

#include "effectivenessHeatExchanger.H"
#include "fvMatrix.H"
#include "basicThermo.H"
#include "surfaceInterpolate.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(effectivenessHeatExchanger, 0);
    addToRunTimeSelectionTable(fvModel, effectivenessHeatExchanger, dictionary);
    addBackwardCompatibleToRunTimeSelectionTable
    (
        fvModel,
        effectivenessHeatExchanger,
        dictionary,
        effectivenessHeatExchangerSource,
        "effectivenessHeatExchangerSource"
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::effectivenessHeatExchanger::readCoeffs()
{
    secondaryMassFlowRate_ = coeffs().lookup<scalar>("secondaryMassFlowRate");
    secondaryInletT_ = coeffs().lookup<scalar>("secondaryInletT");
    primaryInletT_ = coeffs().lookup<scalar>("primaryInletT");

    eTable_.reset(Function2<scalar>::New("effectiveness", coeffs()).ptr());

    UName_ = coeffs().lookupOrDefault<word>("U", "U");
    TName_ = coeffs().lookupOrDefault<word>("T", "T");
    phiName_ = coeffs().lookupOrDefault<word>("phi", "phi");

    faceZoneName_ = coeffs().lookup<word>("faceZone");
}


void Foam::fv::effectivenessHeatExchanger::setZone()
{
    zoneID_ = mesh().faceZones().findZoneID(faceZoneName_);
    if (zoneID_ < 0)
    {
        FatalErrorInFunction
            << type() << " " << this->name() << ": "
            << "    Unknown face zone name: " << faceZoneName_
            << ". Valid face zones are: " << mesh().faceZones().names()
            << nl << exit(FatalError);
    }

    const faceZone& fZone = mesh().faceZones()[zoneID_];

    faceId_.setSize(fZone.size());
    facePatchId_.setSize(fZone.size());
    faceSign_.setSize(fZone.size());

    label count = 0;
    forAll(fZone, i)
    {
        const label facei = fZone[i];
        label faceId = -1;
        label facePatchId = -1;
        if (mesh().isInternalFace(facei))
        {
            faceId = facei;
            facePatchId = -1;
        }
        else
        {
            facePatchId = mesh().boundaryMesh().whichPatch(facei);
            const polyPatch& pp = mesh().boundaryMesh()[facePatchId];
            if (isA<coupledPolyPatch>(pp))
            {
                if (refCast<const coupledPolyPatch>(pp).owner())
                {
                    faceId = pp.whichFace(facei);
                }
                else
                {
                    faceId = -1;
                }
            }
            else if (!isA<emptyPolyPatch>(pp))
            {
                faceId = facei - pp.start();
            }
            else
            {
                faceId = -1;
                facePatchId = -1;
            }
        }

        if (faceId >= 0)
        {
            if (fZone.flipMap()[i])
            {
                faceSign_[count] = -1;
            }
            else
            {
                faceSign_[count] = 1;
            }
            faceId_[count] = faceId;
            facePatchId_[count] = facePatchId;
            count++;
        }
    }
    faceId_.setSize(count);
    facePatchId_.setSize(count);
    faceSign_.setSize(count);

    calculateTotalArea(faceZoneArea_);
}


void Foam::fv::effectivenessHeatExchanger::calculateTotalArea
(
    scalar& area
) const
{
    area = 0;
    forAll(faceId_, i)
    {
        const label facei = faceId_[i];
        if (facePatchId_[i] != -1)
        {
            label patchi = facePatchId_[i];
            area += mesh().magSf().boundaryField()[patchi][facei];
        }
        else
        {
            area += mesh().magSf()[facei];
        }
    }
    reduce(area, sumOp<scalar>());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::effectivenessHeatExchanger::effectivenessHeatExchanger
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvModel(name, modelType, mesh, dict),
    set_(mesh, coeffs()),
    secondaryMassFlowRate_(NaN),
    secondaryInletT_(NaN),
    primaryInletT_(NaN),
    eTable_(nullptr),
    UName_(word::null),
    TName_(word::null),
    phiName_(word::null),
    faceZoneName_(word::null),
    zoneID_(-1),
    faceId_(),
    facePatchId_(),
    faceSign_(),
    faceZoneArea_(NaN)
{
    readCoeffs();
    setZone();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::effectivenessHeatExchanger::addSupFields() const
{
    const basicThermo& thermo =
        mesh().lookupObject<basicThermo>(physicalProperties::typeName);

    return wordList(1, thermo.he().name());
}


void Foam::fv::effectivenessHeatExchanger::addSup
(
    const volScalarField& rho,
    const volScalarField& he,
    fvMatrix<scalar>& eqn
) const
{
    const basicThermo& thermo =
        mesh().lookupObject<basicThermo>(physicalProperties::typeName);

    const surfaceScalarField Cpf(fvc::interpolate(thermo.Cp()));

    const surfaceScalarField& phi =
        mesh().lookupObject<surfaceScalarField>(phiName_);

    scalar totalphi = 0;
    scalar CpfMean = 0;
    forAll(faceId_, i)
    {
        label facei = faceId_[i];
        if (facePatchId_[i] != -1)
        {
            label patchi = facePatchId_[i];
            totalphi += phi.boundaryField()[patchi][facei]*faceSign_[i];

            CpfMean +=
                Cpf.boundaryField()[patchi][facei]
               *mesh().magSf().boundaryField()[patchi][facei];
        }
        else
        {
            totalphi += phi[facei]*faceSign_[i];
            CpfMean += Cpf[facei]*mesh().magSf()[facei];
        }
    }
    reduce(CpfMean, sumOp<scalar>());
    reduce(totalphi, sumOp<scalar>());

    const scalar Qt =
        eTable_->value(mag(totalphi), secondaryMassFlowRate_)
       *(secondaryInletT_ - primaryInletT_)
       *(CpfMean/faceZoneArea_)*mag(totalphi);

    const labelUList cells = set_.cells();

    const volScalarField& T = mesh().lookupObject<volScalarField>(TName_);
    const scalarField TCells(T, cells);
    scalar Tref = 0;
    if (Qt > 0)
    {
        Tref = max(TCells);
        reduce(Tref, maxOp<scalar>());
    }
    else
    {
        Tref = min(TCells);
        reduce(Tref, minOp<scalar>());
    }

    scalarField deltaTCells(cells.size(), 0);
    forAll(deltaTCells, i)
    {
        if (Qt > 0)
        {
            deltaTCells[i] = max(Tref - TCells[i], 0.0);
        }
        else
        {
            deltaTCells[i] = max(TCells[i] - Tref, 0.0);
        }
    }

    const volVectorField& U = mesh().lookupObject<volVectorField>(UName_);
    const scalarField& V = mesh().V();
    scalar sumWeight = 0;

    forAll(cells, i)
    {
        sumWeight += V[cells[i]]*mag(U[cells[i]])*deltaTCells[i];
    }
    reduce(sumWeight, sumOp<scalar>());

    if (mag(Qt) > vSmall)
    {
        scalarField& heSource = eqn.source();

        forAll(cells, i)
        {
            heSource[cells[i]] -=
                Qt*V[cells[i]]*mag(U[cells[i]])*deltaTCells[i]/sumWeight;
        }
    }

    if (debug && Pstream::master())
    {
        Info<< indent << "Net mass flux [Kg/s] = " << totalphi << nl;
        Info<< indent << "Total energy exchange [W] = " << Qt << nl;
        Info<< indent << "Tref [K] = " << Tref << nl;
        Info<< indent << "Effectiveness : "
            << eTable_->value(mag(totalphi), secondaryMassFlowRate_) << endl;
    }
}


bool Foam::fv::effectivenessHeatExchanger::movePoints()
{
    set_.movePoints();
    return true;
}


void Foam::fv::effectivenessHeatExchanger::topoChange
(
    const polyTopoChangeMap& map
)
{
    set_.topoChange(map);
}


void Foam::fv::effectivenessHeatExchanger::mapMesh(const polyMeshMap& map)
{
    set_.mapMesh(map);
}


void Foam::fv::effectivenessHeatExchanger::distribute
(
    const polyDistributionMap& map
)
{
    set_.distribute(map);
}


bool Foam::fv::effectivenessHeatExchanger::read(const dictionary& dict)
{
    if (fvModel::read(dict))
    {
        set_.read(coeffs());
        readCoeffs();
        setZone();
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
