/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "meanVelocityForce.H"
#include "fvMatrices.H"
#include "DimensionedField.H"
#include "IFstream.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(meanVelocityForce, 0);

    addToRunTimeSelectionTable
    (
        option,
        meanVelocityForce,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::meanVelocityForce::readCoeffs()
{
    UName_ = coeffs_.lookupOrDefault<word>("U", "U");

    Ubar_ = coeffs_.lookup<vector>("Ubar");

    relaxation_ = coeffs_.lookupOrDefault<scalar>("relaxation", 1);
}


void Foam::fv::meanVelocityForce::writeProps
(
    const scalar gradP
) const
{
    // Only write on output time
    if (mesh_.time().writeTime())
    {
        IOdictionary propsDict
        (
            IOobject
            (
                name_ + "Properties",
                mesh_.time().timeName(),
                "uniform",
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            )
        );
        propsDict.add("gradient", gradP);
        propsDict.regIOobject::write();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::meanVelocityForce::meanVelocityForce
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    cellSetOption(sourceName, modelType, dict, mesh),
    UName_(word::null),
    Ubar_(vector::uniform(NaN)),
    relaxation_(NaN),
    gradP0_(0),
    dGradP_(0),
    rAPtr_(nullptr)
{
    // Read the initial pressure gradient from file if it exists
    IFstream propsFile
    (
        mesh_.time().timePath()/"uniform"/(name_ + "Properties")
    );
    if (propsFile.good())
    {
        Info<< "    Reading pressure gradient from file" << endl;
        dictionary propsDict(dictionary::null, propsFile);
        propsDict.lookup("gradient") >> gradP0_;
    }

    Info<< "    Initial pressure gradient = " << gradP0_ << nl << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::meanVelocityForce::addSupFields() const
{
    return wordList(1, UName_);
}


Foam::scalar Foam::fv::meanVelocityForce::magUbarAve
(
    const volVectorField& U
) const
{
    const labelList& cells = this->cells();
    const scalarField& cv = mesh_.V();

    scalar magUbarAve = 0;
    forAll(cells, i)
    {
        const label celli = cells[i];
        magUbarAve += (normalised(Ubar_) & U[celli])*cv[celli];
    }
    reduce(magUbarAve, sumOp<scalar>());
    magUbarAve /= V();

    return magUbarAve;
}


void Foam::fv::meanVelocityForce::correct(volVectorField& U) const
{
    const scalarField& rAU = rAPtr_();

    const labelList& cells = this->cells();
    const scalarField& cv = mesh_.V();

    // Average rAU over the cell set
    scalar rAUave = 0;
    forAll(cells, i)
    {
        const label celli = cells[i];
        rAUave += rAU[celli]*cv[celli];
    }
    reduce(rAUave, sumOp<scalar>());
    rAUave /= V();

    const scalar magUbarAve = this->magUbarAve(U);

    // Calculate the pressure gradient increment needed to adjust the average
    // flow-rate to the desired value
    dGradP_ = relaxation_*(mag(Ubar_) - magUbarAve)/rAUave;

    // Apply correction to velocity field
    forAll(cells, i)
    {
        label celli = cells[i];
        U[celli] += normalised(Ubar_)*rAU[celli]*dGradP_;
    }

    const scalar gradP = gradP0_ + dGradP_;

    Info<< "Pressure gradient source: uncorrected Ubar = " << magUbarAve
        << ", pressure gradient = " << gradP << endl;

    writeProps(gradP);
}


void Foam::fv::meanVelocityForce::addSup
(
    fvMatrix<vector>& eqn,
    const word& fieldName
) const
{
    volVectorField::Internal Su
    (
        IOobject
        (
            name_ + fieldName + "Sup",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector(eqn.dimensions()/dimVolume, Zero)
    );

    const scalar gradP = gradP0_ + dGradP_;

    UIndirectList<vector>(Su, cells()) = normalised(Ubar_)*gradP;

    eqn += Su;
}


void Foam::fv::meanVelocityForce::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const word& fieldName
) const
{
    this->addSup(eqn, fieldName);
}


void Foam::fv::meanVelocityForce::constrain
(
    fvMatrix<vector>& eqn,
    const word& fieldName
) const
{
    if (rAPtr_.empty())
    {
        rAPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    name_ + ":rA",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                1/eqn.A()
            )
        );
    }
    else
    {
        rAPtr_() = 1/eqn.A();
    }

    gradP0_ += dGradP_;
    dGradP_ = 0;
}


bool Foam::fv::meanVelocityForce::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        return true;
    }
    else
    {
        return false;
    }
}

// ************************************************************************* //
