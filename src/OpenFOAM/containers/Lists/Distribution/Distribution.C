/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "Distribution.H"
#include "OFstream.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Distribution<Type>::Distribution()
:
    List< List<scalar> >(pTraits<Type>::nComponents),
    binWidth_(pTraits<Type>::one),
    listStarts_(pTraits<Type>::nComponents, 0)
{}


template<class Type>
Foam::Distribution<Type>::Distribution(const Type& binWidth)
:
    List< List<scalar> >(pTraits<Type>::nComponents),
    binWidth_(binWidth),
    listStarts_(pTraits<Type>::nComponents, 0)
{}


template<class Type>
Foam::Distribution<Type>::Distribution(const Distribution<Type>& d)
:
    List< List<scalar> >(static_cast< const List< List<scalar> >& >(d)),
    binWidth_(d.binWidth()),
    listStarts_(d.listStarts())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Distribution<Type>::~Distribution()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
Foam::scalar Foam::Distribution<Type>::totalWeight(direction cmpt) const
{
    const List<scalar>& cmptDistribution = (*this)[cmpt];

    scalar sumOfWeights = 0.0;

    forAll(cmptDistribution, i)
    {
        sumOfWeights += cmptDistribution[i];
    }

    return sumOfWeights;
}


template<class Type>
Foam::List<Foam::label> Foam::Distribution<Type>::keys(direction cmpt) const
{
    List<label> keys = identity((*this)[cmpt].size());

    forAll(keys, k)
    {
        keys[k] += listStarts_[cmpt];
    }

    return keys;
}


template<class Type>
Foam::label Foam::Distribution<Type>::index
(
    direction cmpt,
    label n
)
{
    List<scalar>& cmptDistribution = (*this)[cmpt];

    if (cmptDistribution.empty())
    {
        // Initialise this list with this value

        cmptDistribution.setSize(2, 0.0);

        listStarts_[cmpt] = n;

        return 0;
    }

    label listIndex = -1;

    label& listStart  = listStarts_[cmpt];

    label testIndex = n - listStart;

    if (testIndex < 0)
    {
        // Underflow of this List, storage increase and remapping
        // required

        List<scalar> newCmptDistribution(2*cmptDistribution.size(), 0.0);

        label sOld = cmptDistribution.size();

        forAll(cmptDistribution, i)
        {
            newCmptDistribution[i + sOld] = cmptDistribution[i];
        }

        cmptDistribution = newCmptDistribution;

        listStart -= sOld;

        // Recursively call this function in case another remap is required.
        listIndex = index(cmpt, n);
    }
    else if (testIndex > cmptDistribution.size() - 1)
    {
        // Overflow of this List, storage increase required

        cmptDistribution.setSize(2*cmptDistribution.size(), 0.0);

        // Recursively call this function in case another storage
        // alteration is required.
        listIndex = index(cmpt, n);
    }
    else
    {
        listIndex = n - listStart;
    }

    return listIndex;
}


template<class Type>
Foam::Pair<Foam::label> Foam::Distribution<Type>::validLimits
(
    direction cmpt
) const
{
    const List<scalar>& cmptDistribution = (*this)[cmpt];

    // limits.first(): lower bound, i.e. the first non-zero entry
    // limits.second(): upper bound, i.e. the last non-zero entry
    Pair<label> limits(-1, -1);

    forAll(cmptDistribution, i)
    {
        if (cmptDistribution[i] > 0.0)
        {
            if (limits.first() == -1)
            {
                limits.first() = i;
                limits.second() = i;
            }
            else
            {
                limits.second() = i;
            }
        }
    }

    return limits;
}


template<class Type>
Type Foam::Distribution<Type>::mean() const
{
    Type meanValue(pTraits<Type>::zero);

    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        const List<scalar>& cmptDistribution = (*this)[cmpt];

        scalar totalCmptWeight = totalWeight(cmpt);

        List<label> theKeys = keys(cmpt);

        forAll(theKeys, k)
        {
            label key = theKeys[k];

            setComponent(meanValue, cmpt) +=
                (0.5 + scalar(key))
               *component(binWidth_, cmpt)
               *cmptDistribution[k]
               /totalCmptWeight;
        }
    }

    return meanValue;
}


template<class Type>
Type Foam::Distribution<Type>::median() const
{
    Type medianValue(pTraits<Type>::zero);

    List< List < Pair<scalar> > > normDistribution = normalised();

    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        List< Pair<scalar> >& normDist = normDistribution[cmpt];

        if (normDist.size())
        {
            if (normDist.size() == 1)
            {
                setComponent(medianValue, cmpt) = normDist[0].first();
            }
            else if
            (
                normDist.size() > 1
             && normDist[0].second()*component(binWidth_, cmpt) > 0.5
            )
            {
                scalar xk =
                    normDist[0].first()
                  + 0.5*component(binWidth_, cmpt);

                scalar xkm1 =
                    normDist[0].first()
                  - 0.5*component(binWidth_, cmpt);

                scalar Sk = (normDist[0].second())*component(binWidth_, cmpt);

                setComponent(medianValue, cmpt) = 0.5*(xk - xkm1)/(Sk) + xkm1;
            }
            else
            {
                label previousNonZeroIndex = 0;

                scalar cumulative = 0.0;

                forAll(normDist, nD)
                {
                    if
                    (
                        cumulative
                      + (normDist[nD].second()*component(binWidth_, cmpt))
                      > 0.5
                    )
                    {
                        scalar xk =
                            normDist[nD].first()
                          + 0.5*component(binWidth_, cmpt);

                        scalar xkm1 =
                            normDist[previousNonZeroIndex].first()
                          + 0.5*component(binWidth_, cmpt);

                        scalar Sk =
                            cumulative
                          + (normDist[nD].second()*component(binWidth_, cmpt));

                        scalar Skm1 = cumulative;

                        setComponent(medianValue, cmpt) =
                            (0.5 - Skm1)*(xk - xkm1)/(Sk - Skm1) + xkm1;

                        break;
                    }
                    else if (mag(normDist[nD].second()) > VSMALL)
                    {
                        cumulative +=
                            normDist[nD].second()*component(binWidth_, cmpt);

                        previousNonZeroIndex = nD;
                    }
                }
            }
        }

    }

    return medianValue;
}


template<class Type>
void Foam::Distribution<Type>::add
(
    const Type& valueToAdd,
    const Type& weight
)
{
    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        List<scalar>& cmptDistribution = (*this)[cmpt];

        label n =
            label(component(valueToAdd, cmpt)/component(binWidth_, cmpt))
          - label(neg(component(valueToAdd, cmpt)/component(binWidth_, cmpt)));

        label listIndex = index(cmpt, n);

        cmptDistribution[listIndex] += component(weight, cmpt);
    }
}


template<class Type>
Foam::List< Foam::List< Foam::Pair<Foam::scalar> > >Foam::
Distribution<Type>::normalised() const
{
    List< List < Pair<scalar> > > normDistribution(pTraits<Type>::nComponents);

    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        const List<scalar>& cmptDistribution = (*this)[cmpt];

        if (cmptDistribution.empty())
        {
            continue;
        }

        scalar totalCmptWeight = totalWeight(cmpt);

        List<label> cmptKeys = keys(cmpt);

        List< Pair<scalar> >& normDist = normDistribution[cmpt];

        Pair<label> limits = validLimits(cmpt);

        normDist.setSize(limits.second() - limits.first() + 1);

        for
        (
            label k = limits.first(), i = 0;
            k <= limits.second();
            k++, i++
        )
        {
            label key = cmptKeys[k];

            normDist[i].first() =
                (0.5 + scalar(key))*component(binWidth_, cmpt);

            normDist[i].second() =
                cmptDistribution[k]
               /totalCmptWeight
               /component(binWidth_, cmpt);
        }
    }

    return normDistribution;
}


template<class Type>
Foam::List< Foam::List< Foam::Pair<Foam::scalar> > >Foam::
Distribution<Type>::raw() const
{
    List< List < Pair<scalar> > > rawDistribution(pTraits<Type>::nComponents);

    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        const List<scalar>& cmptDistribution = (*this)[cmpt];

        if (cmptDistribution.empty())
        {
            continue;
        }

        List<label> cmptKeys = keys(cmpt);

        List< Pair<scalar> >& rawDist = rawDistribution[cmpt];

        Pair<label> limits = validLimits(cmpt);

        rawDist.setSize(limits.second() - limits.first() + 1);

        for
        (
            label k = limits.first(), i = 0;
            k <= limits.second();
            k++, i++
        )
        {
            label key = cmptKeys[k];

            rawDist[i].first() = (0.5 + scalar(key))*component(binWidth_, cmpt);

            rawDist[i].second() = cmptDistribution[k];
        }
    }

    return rawDistribution;
}


template<class Type>
Foam::List< Foam::List< Foam::Pair<Foam::scalar> > >Foam::
Distribution<Type>::cumulativeNormalised() const
{
    List< List< Pair<scalar> > > normalisedDistribution = normalised();

    List< List < Pair<scalar> > > cumulativeNormalisedDistribution =
        normalisedDistribution;

    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        const List< Pair<scalar> >& normalisedCmpt =
            normalisedDistribution[cmpt];

        List< Pair<scalar> >& cumNormalisedCmpt =
            cumulativeNormalisedDistribution[cmpt];

        scalar sum = 0.0;

        forAll(normalisedCmpt, i)
        {
            cumNormalisedCmpt[i].first() =
                normalisedCmpt[i].first()
              + 0.5*component(binWidth_, cmpt);

            cumNormalisedCmpt[i].second() =
                normalisedCmpt[i].second()*component(binWidth_, cmpt) + sum;

            sum = cumNormalisedCmpt[i].second();
        }
    }

    return cumulativeNormalisedDistribution;
}


template<class Type>
Foam::List< Foam::List< Foam::Pair<Foam::scalar> > >Foam::
Distribution<Type>::cumulativeRaw() const
{
    List< List< Pair<scalar> > > rawDistribution = raw();

    List< List < Pair<scalar> > > cumulativeRawDistribution = rawDistribution;

    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        const List< Pair<scalar> >& rawCmpt = rawDistribution[cmpt];

        List< Pair<scalar> >& cumRawCmpt = cumulativeRawDistribution[cmpt];

        scalar sum = 0.0;

        forAll(rawCmpt, i)
        {
            cumRawCmpt[i].first() =
                rawCmpt[i].first()
              + 0.5*component(binWidth_, cmpt);

            cumRawCmpt[i].second() = rawCmpt[i].second() + sum;

            sum = cumRawCmpt[i].second();
        }
    }

    return cumulativeRawDistribution;
}


template<class Type>
void Foam::Distribution<Type>::clear()
{
    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        (*this)[cmpt].clear();

        listStarts_[cmpt] = 0;
    }
}


template<class Type>
void Foam::Distribution<Type>::write(const fileName& filePrefix) const
{
    List< List< Pair<scalar> > > rawDistribution = raw();

    List< List < Pair<scalar> > > normDistribution = normalised();

    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        const List< Pair<scalar> >& rawPairs = rawDistribution[cmpt];

        const List< Pair<scalar> >& normPairs = normDistribution[cmpt];

        OFstream os(filePrefix + '_' + pTraits<Type>::componentNames[cmpt]);

        os  << "# key normalised raw" << endl;

        forAll(normPairs, i)
        {
            os  << normPairs[i].first()
                << ' ' << normPairs[i].second()
                << ' ' << rawPairs[i].second()
                << nl;
        }
    }

    List< List< Pair<scalar> > > rawCumDist = cumulativeRaw();

    List< List < Pair<scalar> > > normCumDist = cumulativeNormalised();

    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        const List< Pair<scalar> >& rawPairs = rawCumDist[cmpt];

        const List< Pair<scalar> >& normPairs = normCumDist[cmpt];

        OFstream os
        (
            filePrefix + "_cumulative_" + pTraits<Type>::componentNames[cmpt]
        );

        os  << "# key normalised raw" << endl;

        forAll(normPairs, i)
        {
            os  << normPairs[i].first()
                << ' ' << normPairs[i].second()
                << ' ' << rawPairs[i].second()
                << nl;
        }
    }
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

template<class Type>
void Foam::Distribution<Type>::operator=
(
    const Distribution<Type>& rhs
)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn
        (
            "Foam::Distribution<Type>::operator="
            "(const Foam::Distribution<Type>&)"
        )   << "Attempted assignment to self"
            << abort(FatalError);
    }

    List< List<scalar> >::operator=(rhs);

    binWidth_ = rhs.binWidth();

    listStarts_ = rhs.listStarts();
}


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

template<class Type>
Foam::Istream& Foam::operator>>
(
    Istream& is,
    Distribution<Type>& d
)
{
    is  >> static_cast<List< List<scalar> >&>(d)
        >> d.binWidth_
        >> d.listStarts_;

    // Check state of Istream
    is.check("Istream& operator>>(Istream&, Distribution<Type>&)");

    return is;
}


template<class Type>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const Distribution<Type>& d
)
{
    os  <<  static_cast<const List< List<scalar> >& >(d)
        << d.binWidth_ << token::SPACE
        << d.listStarts_;

    // Check state of Ostream
    os.check("Ostream& operator<<(Ostream&, " "const Distribution&)");

    return os;
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

template<class Type>
Foam::Distribution<Type> Foam::operator+
(
    const Distribution<Type>& d1,
    const Distribution<Type>& d2
)
{
    // The coarsest binWidth is the sensible choice
    Distribution<Type> d(max(d1.binWidth(), d2.binWidth()));

    List< List< List < Pair<scalar> > > > rawDists(2);

    rawDists[0] = d1.raw();
    rawDists[1] = d2.raw();

    forAll(rawDists, rDI)
    {
        for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
        {
            List<scalar>& cmptDistribution = d[cmpt];

            const List < Pair<scalar> >& cmptRaw = rawDists[rDI][cmpt];

            forAll(cmptRaw, rI)
            {
                scalar valueToAdd = cmptRaw[rI].first();
                scalar cmptWeight = cmptRaw[rI].second();

                label n =
                label
                (
                    component(valueToAdd, cmpt)
                   /component(d.binWidth(), cmpt)
                )
                - label
                (
                    neg(component(valueToAdd, cmpt)
                   /component(d.binWidth(), cmpt))
                );

                label listIndex = d.index(cmpt, n);

                cmptDistribution[listIndex] += cmptWeight;
            }
        }
    }

    return Distribution<Type>(d);
}


// ************************************************************************* //
