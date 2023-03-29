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

#include "tabulatedCumulative.H"
#include "unintegrable.H"
#include "SubField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace distributions
{
    defineTypeNameAndDebug(tabulatedCumulative, 0);
    addToRunTimeSelectionTable(distribution, tabulatedCumulative, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::distributions::tabulatedCumulative::tabulatedCumulative
(
    const dictionary& dict,
    Random& rndGen,
    const label sampleQ
)
:
    FieldDistribution<distribution, tabulatedCumulative>
    (
        typeName,
        dict,
        rndGen,
        sampleQ
    )
{
    List<Pair<scalar>> values(dict.lookup("distribution"));

    // Checks
    if (values.first().second() != 0)
    {
        FatalIOErrorInFunction(dict)
            << typeName << ": The cumulative distribution function does not "
            << "start at zero" << abort(FatalIOError);
    }

    for (label i = 1; i < values.size(); ++ i)
    {
        if (values[i].first() - values[i - 1].first() < 0)
        {
            FatalIOErrorInFunction(dict)
                << typeName << ": The cumulative distribution function is not "
                << "in order" << abort(FatalIOError);
        }

        if (values[i].second() - values[i - 1].second() < 0)
        {
            FatalIOErrorInFunction(dict)
                << typeName << ": The cumulative distribution function is not "
                << "monotonic" << abort(FatalIOError);
        }
    }

    // Copy the coordinates
    x_.resize(values.size());
    forAll(values, i)
    {
        x_[i] = values[i].first();
    }

    // Set the CDF. Copy if q == 0. Re-integrated if q != 0.
    CDF_.resize(values.size());
    CDF_[0] = 0;
    for (label i = 1; i < values.size(); ++ i)
    {
        CDF_[i] =
            CDF_[i - 1]
          + integerPow((x_[i] + x_[i - 1])/2, q())
           *(values[i].second() - values[i - 1].second());
    }
    CDF_ /= CDF_.last();

    // Compute the PDF. Equals the gradient of the CDF. Associated with the
    // intervals, so the field is one element shorter than the coordinates and
    // the CDF.
    PDF_.resize(values.size() - 1);
    forAll(PDF_, i)
    {
        PDF_[i] = (CDF_[i + 1] - CDF_[i])/(x_[i + 1] - x_[i]);
    }

    // More checks
    validateBounds(dict);
    if (q() != 0) validatePositive(dict);
    report();
}


Foam::distributions::tabulatedCumulative::tabulatedCumulative
(
    const tabulatedCumulative& d,
    const label sampleQ
)
:
    FieldDistribution<distribution, tabulatedCumulative>(d, sampleQ),
    x_(d.x_),
    PDF_(d.PDF_),
    CDF_(d.CDF_)
{
    // If Q is the same then a copy is sufficient
    if (q() == d.q()) return;

    // Re-integrate the CDF
    CDF_[0] = 0;
    for (label i = 1; i < x_.size(); ++ i)
    {
        CDF_[i] =
            CDF_[i - 1]
          + integerPow((d.x_[i] + d.x_[i - 1])/2, q() - d.q())
           *(d.CDF_[i] - d.CDF_[i - 1]);
    }
    CDF_ /= CDF_.last();

    // Re-compute the PDF
    forAll(PDF_, i)
    {
        PDF_[i] = (CDF_[i + 1] - CDF_[i])/(x_[i + 1] - x_[i]);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::distributions::tabulatedCumulative::~tabulatedCumulative()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::distributions::tabulatedCumulative::sample() const
{
    return unintegrable::sample(x_, CDF_, rndGen_.sample01<scalar>());
}


Foam::scalar Foam::distributions::tabulatedCumulative::min() const
{
    return x_.first();
}


Foam::scalar Foam::distributions::tabulatedCumulative::max() const
{
    return x_.last();
}


Foam::scalar Foam::distributions::tabulatedCumulative::mean() const
{
    const scalarField x(this->x(-1));
    return unintegrable::integrate(x, x*PDF(x))->last();
}


Foam::tmp<Foam::scalarField> Foam::distributions::tabulatedCumulative::x
(
    const label
) const
{
    const scalar x0 = min(), x1 = max(), d = 0.1*(x1 - x0);

    tmp<scalarField> tResult(new scalarField(2*x_.size() + 2));
    scalarField& result = tResult.ref();

    result[0] = Foam::max(x0 - d, q() < 0 ? x0/2 : rootVSmall);

    forAll(x_, i)
    {
        result[2*i + 1] = result[2*i + 2] = x_[i];
    }

    result[2*x_.size() + 1] = x1 + d;

    return tResult;
}


Foam::tmp<Foam::scalarField>
Foam::distributions::tabulatedCumulative::PDF(const scalarField& x) const
{
    tmp<scalarField> tResult(new scalarField(2*PDF_.size() + 4));
    scalarField& result = tResult.ref();

    result[0] = result[1] = 0;

    forAll(PDF_, i)
    {
        result[2*i + 2] = result[2*i + 3] = PDF_[i];
    }

    result[2*PDF_.size() + 2] = result[2*PDF_.size() + 3] = 0;

    return tResult;
}


// ************************************************************************* //
