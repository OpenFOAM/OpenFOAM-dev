/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "distribution.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(distribution, 0);
    defineRunTimeSelectionTable(distribution, dictionary);
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::distribution::validateBounds(const dictionary& dict) const
{
    if (max() < min())
    {
        FatalIOErrorInFunction(dict)
            << type() << ": The maximum value is smaller than the minimum "
            << "value:" << nl << "    max = " << max() << ", min = "
            << min() << abort(FatalIOError);
    }
}


void Foam::distribution::validatePositive(const dictionary& dict) const
{
    if (min() < 0)
    {
        FatalIOErrorInFunction(dict)
            << type() << ": The minimum value must be greater than "
            << "zero." << nl << "    min = " << min()
            << abort(FatalIOError);
    }
}


Foam::tmp<Foam::scalarField>
Foam::distribution::clipPDF
(
    const scalarField& x,
    const tmp<scalarField>& pdf
) const
{
    return pos0(x - min())*pos0(max() - x)*pdf;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::distribution::distribution
(
    const word& name,
    const unitConversion& units,
    const dictionary& dict,
    const label sampleQ,
    randomGenerator&& rndGen
)
:
    Q_(dict.lookup<scalar>("Q")),
    sampleQ_(sampleQ),
    rndGen_("rndGen", dict, std::move(rndGen))
{
    if (Q_ < 0)
    {
        FatalIOErrorInFunction(dict)
            << name << ": Size exponent cannot be negative" << nl
            << "    Q = " << Q_ << abort(FatalIOError);
    }

    if (sampleQ_ < 0)
    {
        FatalErrorInFunction
            << name << ": Sampling size exponent cannot be negative" << nl
            << "    sampleQ = " << sampleQ_ << abort(FatalError);
    }
}


Foam::distribution::distribution
(
    const label Q,
    const label sampleQ,
    randomGenerator&& rndGen
)
:
    Q_(Q),
    sampleQ_(sampleQ),
    rndGen_(rndGen)
{}


Foam::distribution::distribution(const distribution& d, const label sampleQ)
:
    Q_(d.Q_),
    sampleQ_(sampleQ),
    rndGen_(d.rndGen_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::distribution::~distribution()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::distribution::write(Ostream& os, const unitConversion& units) const
{
    writeEntry(os, "type", type());
    writeEntry(os, "Q", Q_);
}


void Foam::distribution::writeState(Ostream& os) const
{
    writeEntry(os, "rndGen", rndGen_);
}


Foam::tmp<Foam::scalarField> Foam::distribution::plotX(const label n) const
{
    const scalar x0 = min(), x1 = max(), d = 0.1*(x1 - x0);

    tmp<scalarField> tResult(new scalarField(n));
    scalarField& result = tResult.ref();

    result[0] = x0 - d;
    result[1] = x0*(1 - sign(x0)*small);
    result[2] = x0*(1 + sign(x0)*small);

    for (label i = 3; i < n - 3; ++ i)
    {
        const scalar f = scalar(i - 2)/(n - 5);

        result[i] = (1 - f)*x0 + f*x1;
    }

    result[n - 3] = x1*(1 - sign(x1)*small);
    result[n - 2] = x1*(1 + sign(x1)*small);
    result[n - 1] = x1 + d;

    return tResult;
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void Foam::writeEntry
(
    Ostream& os,
    const word& entryName,
    const unitConversion& units,
    const distribution& d,
    const bool write,
    const bool writeState
)
{
    writeKeyword(os, entryName);

    os  << nl << indent << token::BEGIN_BLOCK << nl << incrIndent;

    if (write) d.write(os, units);
    if (writeState) d.writeState(os);

    os  << decrIndent << indent << token::END_BLOCK << endl;
}


// ************************************************************************* //
