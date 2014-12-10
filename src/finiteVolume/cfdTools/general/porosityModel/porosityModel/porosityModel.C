/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2014 OpenFOAM Foundation
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

#include "porosityModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(porosityModel, 0);
    defineRunTimeSelectionTable(porosityModel, mesh);
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::porosityModel::adjustNegativeResistance(dimensionedVector& resist)
{
    scalar maxCmpt = max(0, cmptMax(resist.value()));

    if (maxCmpt < 0)
    {
        FatalErrorIn
        (
            "void Foam::porosityModel::adjustNegativeResistance"
            "("
                "dimensionedVector&"
            ")"
        )   << "Negative resistances are invalid, resistance = " << resist
            << exit(FatalError);
    }
    else
    {
        vector& val = resist.value();
        for (label cmpt = 0; cmpt < vector::nComponents; cmpt++)
        {
            if (val[cmpt] < 0)
            {
                val[cmpt] *= -maxCmpt;
            }
        }
    }
}


Foam::label Foam::porosityModel::fieldIndex(const label i) const
{
    label index = 0;
    if (!coordSys_.R().uniform())
    {
        index = i;
    }
    return index;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porosityModel::porosityModel
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& cellZoneName
)
:
    regIOobject
    (
        IOobject
        (
            name,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    name_(name),
    mesh_(mesh),
    dict_(dict),
    coeffs_(dict.subDict(modelType + "Coeffs")),
    active_(true),
    zoneName_(cellZoneName),
    cellZoneIDs_(),
    coordSys_(coordinateSystem::New(mesh, coeffs_))
{
    if (zoneName_ == word::null)
    {
        dict.lookup("active") >> active_;
        dict_.lookup("cellZone") >> zoneName_;
    }

    cellZoneIDs_ = mesh_.cellZones().findIndices(zoneName_);

    Info<< "    creating porous zone: " << zoneName_ << endl;

    bool foundZone = !cellZoneIDs_.empty();
    reduce(foundZone, orOp<bool>());

    if (!foundZone && Pstream::master())
    {
        FatalErrorIn
        (
            "Foam::porosityModel::porosityModel"
            "("
                "const word&, "
                "const word&, "
                "const fvMesh&, "
                "const dictionary&"
                "const word&, "
            ")"
        )   << "cannot find porous cellZone " << zoneName_
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::porosityModel::~porosityModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::porosityModel::transformModelData()
{
    if (!mesh_.upToDatePoints(*this))
    {
        calcTranformModelData();

        // set model up-to-date wrt points
        mesh_.setUpToDatePoints(*this);
    }
}


Foam::tmp<Foam::vectorField> Foam::porosityModel::porosityModel::force
(
    const volVectorField& U,
    const volScalarField& rho,
    const volScalarField& mu
)
{
    transformModelData();

    tmp<vectorField> tforce(new vectorField(U.size(), vector::zero));

    if (!cellZoneIDs_.empty())
    {
        this->calcForce(U, rho, mu, tforce());
    }

    return tforce;
}


void Foam::porosityModel::addResistance(fvVectorMatrix& UEqn)
{
    if (cellZoneIDs_.empty())
    {
        return;
    }

    transformModelData();
    this->correct(UEqn);
}


void Foam::porosityModel::addResistance
(
    fvVectorMatrix& UEqn,
    const volScalarField& rho,
    const volScalarField& mu
)
{
    if (cellZoneIDs_.empty())
    {
        return;
    }

    transformModelData();
    this->correct(UEqn, rho, mu);
}


void Foam::porosityModel::addResistance
(
    const fvVectorMatrix& UEqn,
    volTensorField& AU,
    bool correctAUprocBC
)
{
    if (cellZoneIDs_.empty())
    {
        return;
    }

    transformModelData();
    this->correct(UEqn, AU);

    if (correctAUprocBC)
    {
        // Correct the boundary conditions of the tensorial diagonal to ensure
        // processor boundaries are correctly handled when AU^-1 is interpolated
        // for the pressure equation.
        AU.correctBoundaryConditions();
    }
}


bool Foam::porosityModel::writeData(Ostream& os) const
{
    return true;
}


bool Foam::porosityModel::read(const dictionary& dict)
{
    active_ = readBool(dict.lookup("active"));
    coeffs_ = dict.subDict(type() + "Coeffs");
    dict.lookup("cellZone") >> zoneName_;
    cellZoneIDs_ = mesh_.cellZones().findIndices(zoneName_);

    return true;
}


// ************************************************************************* //
