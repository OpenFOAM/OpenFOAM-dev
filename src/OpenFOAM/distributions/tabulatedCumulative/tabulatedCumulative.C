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
    const unitConversion& defaultUnits,
    const dictionary& dict,
    const label sampleQ,
    randomGenerator&& rndGen
)
:
    FieldDistribution<distribution, tabulatedCumulative>
    (
        typeName,
        defaultUnits,
        dict,
        sampleQ,
        std::move(rndGen)
    ),
    reader_
    (
        TableReader<scalar>::New(word::null, {defaultUnits, unitAny}, dict)
    )
{
    List<Tuple2<scalar, scalar>> values =
        reader_->read({defaultUnits, unitAny}, dict, "distribution");

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
}


Foam::distributions::tabulatedCumulative::tabulatedCumulative
(
    const tabulatedCumulative& d,
    const label sampleQ
)
:
    FieldDistribution<distribution, tabulatedCumulative>(d, sampleQ),
    reader_(d.reader_, false),
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
    const scalarField x(this->plotX(-1));
    return unintegrable::integrate(x, x*plotPDF(x))->last();
}


Foam::tmp<Foam::scalarField>
Foam::distributions::tabulatedCumulative::integralPDFxPow
(
    const scalarField& x,
    const label e,
    const bool
) const
{
    const scalarField& xStar = x_;
    const scalarField& yStar = PDF_;

    tmp<scalarField> tResult(new scalarField(x.size()));
    scalarField& result = tResult.ref();

    label i = 0;

    while (i < x.size() && x[i] < xStar[0])
    {
        result[i] = 0;
        i ++;
    }

    scalar integral_PDFxPowE_0_j = 0;

    for (label iStar = 0; iStar < xStar.size() - 1; ++ iStar)
    {
        const scalar xPowE1_j = integerPow(xStar[iStar], e + 1);

        auto integral_xPowE1_j_x = [&](const scalar x)
        {
            const scalar xPowE1_i = integerPow(x, e + 1);

            const scalar integral_xPowE_j_x =
                e + 1 == 0
              ? log(x/xStar[iStar])
              : (xPowE1_i - xPowE1_j)/(e + 1);

            return yStar[iStar]*integral_xPowE_j_x;
        };

        while (i < x.size() && x[i] < xStar[iStar + 1])
        {
            result[i] = integral_PDFxPowE_0_j + integral_xPowE1_j_x(x[i]);

            i ++;
        }

        integral_PDFxPowE_0_j += integral_xPowE1_j_x(xStar[iStar + 1]);
    }

    while (i < x.size())
    {
        result[i] = integral_PDFxPowE_0_j;
        i ++;
    }

    return tResult;

}


void Foam::distributions::tabulatedCumulative::write
(
    Ostream& os,
    const unitConversion& units
) const
{
    FieldDistribution<distribution, tabulatedCumulative>::write(os, units);

    // Recover the CDF
    scalarField CDF(CDF_.size());
    CDF[0] = 0;
    for (label i = 1; i < CDF_.size(); ++ i)
    {
        CDF[i] =
            CDF[i - 1]
          + integerPow((x_[i] + x_[i - 1])/2, -q())
           *(CDF_[i] - CDF_[i - 1]);
    }

    // Normalise the CDF
    CDF /= CDF.last();

    // Construct and write the values
    List<Tuple2<scalar, scalar>> values(CDF_.size());
    forAll(values, i)
    {
        values[i].first() = x_[i];
        values[i].second() = CDF[i];
    }
    reader_->write(os, {units, unitAny}, values, "distribution");
}


Foam::tmp<Foam::scalarField>
Foam::distributions::tabulatedCumulative::plotX(const label) const
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
Foam::distributions::tabulatedCumulative::plotPDF(const scalarField& x) const
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
