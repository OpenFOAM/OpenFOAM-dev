/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2021 OpenFOAM Foundation
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

#include "effectivenessHeatExchangerSource.H"
#include "fvMatrix.H"
#include "basicThermo.H"
#include "surfaceInterpolate.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(effectivenessHeatExchangerSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        effectivenessHeatExchangerSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::effectivenessHeatExchangerSource::readCoeffs()
{
    secondaryMassFlowRate_ = coeffs_.lookup<scalar>("secondaryMassFlowRate");
    secondaryInletT_ = coeffs_.lookup<scalar>("secondaryInletT");
    primaryInletT_ = coeffs_.lookup<scalar>("primaryInletT");

    eTable_.reset(Function2<scalar>::New("effectiveness", coeffs_).ptr());

    UName_ = coeffs_.lookupOrDefault<word>("U", "U");
    TName_ = coeffs_.lookupOrDefault<word>("T", "T");
    phiName_ = coeffs_.lookupOrDefault<word>("phi", "phi");

    faceZoneName_ = coeffs_.lookup<word>("faceZone");
}


void Foam::fv::effectivenessHeatExchangerSource::setZone()
{
    zoneID_ = mesh_.faceZones().findZoneID(faceZoneName_);
    if (zoneID_ < 0)
    {
        FatalErrorInFunction
            << type() << " " << this->name() << ": "
            << "    Unknown face zone name: " << faceZoneName_
            << ". Valid face zones are: " << mesh_.faceZones().names()
            << nl << exit(FatalError);
    }

    const faceZone& fZone = mesh_.faceZones()[zoneID_];

    faceId_.setSize(fZone.size());
    facePatchId_.setSize(fZone.size());
    faceSign_.setSize(fZone.size());

    label count = 0;
    forAll(fZone, i)
    {
        const label facei = fZone[i];
        label faceId = -1;
        label facePatchId = -1;
        if (mesh_.isInternalFace(facei))
        {
            faceId = facei;
            facePatchId = -1;
        }
        else
        {
            facePatchId = mesh_.boundaryMesh().whichPatch(facei);
            const polyPatch& pp = mesh_.boundaryMesh()[facePatchId];
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


void Foam::fv::effectivenessHeatExchangerSource::calculateTotalArea
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
            area += mesh_.magSf().boundaryField()[patchi][facei];
        }
        else
        {
            area += mesh_.magSf()[facei];
        }
    }
    reduce(area, sumOp<scalar>());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::effectivenessHeatExchangerSource::effectivenessHeatExchangerSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    cellSetOption(name, modelType, dict, mesh),
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

Foam::wordList Foam::fv::effectivenessHeatExchangerSource::addSupFields() const
{
    const basicThermo& thermo =
        mesh_.lookupObject<basicThermo>(basicThermo::dictName);

    return wordList(1, thermo.he().name());
}


void Foam::fv::effectivenessHeatExchangerSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    const basicThermo& thermo =
        mesh_.lookupObject<basicThermo>(basicThermo::dictName);

    const surfaceScalarField Cpf(fvc::interpolate(thermo.Cp()));

    const surfaceScalarField& phi =
        mesh_.lookupObject<surfaceScalarField>(phiName_);

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
               *mesh_.magSf().boundaryField()[patchi][facei];
        }
        else
        {
            totalphi += phi[facei]*faceSign_[i];
            CpfMean += Cpf[facei]*mesh_.magSf()[facei];
        }
    }
    reduce(CpfMean, sumOp<scalar>());
    reduce(totalphi, sumOp<scalar>());

    const scalar Qt =
        eTable_->value(mag(totalphi), secondaryMassFlowRate_)
       *(secondaryInletT_ - primaryInletT_)
       *(CpfMean/faceZoneArea_)*mag(totalphi);

    const volScalarField& T = mesh_.lookupObject<volScalarField>(TName_);
    const scalarField TCells(T, cells());
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

    scalarField deltaTCells(cells().size(), 0);
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

    const volVectorField& U = mesh_.lookupObject<volVectorField>(UName_);
    const scalarField& V = mesh_.V();
    scalar sumWeight = 0;

    const labelList& cells = this->cells();

    forAll(cells, i)
    {
        sumWeight += V[cells[i]]*mag(U[cells[i]])*deltaTCells[i];
    }
    reduce(sumWeight, sumOp<scalar>());

    if (this->V() > vSmall && mag(Qt) > vSmall)
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


bool Foam::fv::effectivenessHeatExchangerSource::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
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
