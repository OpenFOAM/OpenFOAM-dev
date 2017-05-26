/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenFOAM Foundation
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

#include "waveModel.H"
#include "mathematicalConstants.H"
#include "Time.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(waveModel, 0);
    defineRunTimeSelectionTable(waveModel, objectRegistry);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::waveModel::k() const
{
    return 2*Foam::constant::mathematical::pi/length_;
}


Foam::scalar Foam::waveModel::sigma() const
{
    const uniformDimensionedVectorField& g =
        db_.lookupObject<uniformDimensionedVectorField>("g");

    return sqrt(mag(g.value())*k()*tanh(k()*depth()));
}


Foam::scalar Foam::waveModel::omega(const scalar u) const
{
    return sigma() + k()*u;
}


Foam::tmp<Foam::scalarField> Foam::waveModel::angle
(
    const scalar t,
    const scalar u,
    const scalarField& x
) const
{
    return k()*x - omega(u)*t;
}


bool Foam::waveModel::shallow() const
{
    return k()*depth() < log(GREAT);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waveModel::waveModel(const waveModel& wave)
:
    db_(wave.db_),
    length_(wave.length_),
    amplitude_(wave.amplitude_, false),
    phase_(wave.phase_),
    depth_(wave.depth_)
{}


Foam::waveModel::waveModel(const objectRegistry& db, const dictionary& dict)
:
    db_(db),
    length_(readScalar(dict.lookup("length"))),
    amplitude_(Function1<scalar>::New("amplitude", dict)),
    phase_(readScalar(dict.lookup("phase"))),
    depth_(dict.lookupOrDefault<scalar>("depth", GREAT))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::waveModel::~waveModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::waveModel::write(Ostream& os) const
{
    os.writeKeyword("length") << length_ << token::END_STATEMENT << nl;
    amplitude_->writeData(os);
    os.writeKeyword("phase") << phase_ << token::END_STATEMENT << nl;
    os.writeKeyword("depth") << depth_ << token::END_STATEMENT << nl;
}


// ************************************************************************* //
