/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "tabulatedDensity.H"
#include "unintegrable.H"
#include "SubField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace distributions
{
    defineTypeNameAndDebug(tabulatedDensity, 0);
    addToRunTimeSelectionTable(distribution, tabulatedDensity, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::distributions::tabulatedDensity::tabulatedDensity
(
    const dictionary& dict,
    Random& rndGen,
    const label sampleQ
)
:
    FieldDistribution<distribution, tabulatedDensity>
    (
        typeName,
        dict,
        rndGen,
        sampleQ
    )
{
    List<Pair<scalar>> values(dict.lookup("distribution"));

    // Checks
    forAll(values, i)
    {
        if (i && values[i].first() - values[i - 1].first() < 0)
        {
            FatalIOErrorInFunction(dict)
                << typeName << ": The probability density function is not "
                << "in order" << abort(FatalIOError);
        }

        if (values[i].second() < 0)
        {
            FatalIOErrorInFunction(dict)
                << typeName << ": The probability density function is not "
                << "positive" << abort(FatalIOError);
        }
    }

    // Copy the coordinates
    x_.resize(values.size());
    forAll(values, i)
    {
        x_[i] = values[i].first();
    }

    // Copy the PDF. Scale if q != 0.
    PDF_.resize(values.size());
    forAll(values, i)
    {
        PDF_[i] = integerPow(x_[i], q())*values[i].second();
    }

    // Compute the CDF
    CDF_.resize(values.size());
    CDF_ = unintegrable::integrate(x_, PDF_);

    // Normalise the PDF and the CDF
    PDF_ /= CDF_.last();
    CDF_ /= CDF_.last();

    // More checks
    validateBounds(dict);
    if (q() != 0) validatePositive(dict);
    report();
}


Foam::distributions::tabulatedDensity::tabulatedDensity
(
    const tabulatedDensity& d,
    const label sampleQ
)
:
    FieldDistribution<distribution, tabulatedDensity>(d, sampleQ),
    x_(d.x_),
    PDF_(d.PDF_),
    CDF_(d.CDF_)
{
    if (q() == d.q()) return;

    // Scale the PDF
    PDF_ = integerPow(x_, q() - d.q())*d.PDF_;

    // Compute the CDF
    CDF_ = unintegrable::integrate(x_, PDF_);

    // Normalise the PDF and the CDF
    PDF_ /= CDF_.last();
    CDF_ /= CDF_.last();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::distributions::tabulatedDensity::~tabulatedDensity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::distributions::tabulatedDensity::sample() const
{
    return unintegrable::sample(x_, PDF_, CDF_, rndGen_.sample01<scalar>());
}


Foam::scalar Foam::distributions::tabulatedDensity::min() const
{
    return x_.first();
}


Foam::scalar Foam::distributions::tabulatedDensity::max() const
{
    return x_.last();
}


Foam::scalar Foam::distributions::tabulatedDensity::mean() const
{
    return unintegrable::integrateX(x_, PDF_)->last();
}


Foam::tmp<Foam::scalarField> Foam::distributions::tabulatedDensity::x
(
    const label
) const
{
    const scalar x0 = min(), x1 = max(), d = 0.1*(x1 - x0);

    tmp<scalarField> tResult(new scalarField(x_.size() + 4));
    scalarField& result = tResult.ref();

    result[0] = Foam::max(x0 - d, q() < 0 ? x0/2 : rootVSmall);
    result[1] = x0;

    SubField<scalar>(result, x_.size(), 2) = x_;

    result[x_.size() + 2] = x1;
    result[x_.size() + 3] = x1 + d;

    return tResult;
}


Foam::tmp<Foam::scalarField>
Foam::distributions::tabulatedDensity::PDF(const scalarField& x) const
{
    tmp<scalarField> tResult(new scalarField(PDF_.size() + 4, 0));
    scalarField& result = tResult.ref();

    SubField<scalar>(result, PDF_.size(), 2) = PDF_;

    return tResult;
}


// ************************************************************************* //
