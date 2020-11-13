/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "interpolationLookUpTable.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::label
Foam::radiationModels::absorptionEmissionModels::interpolationLookUpTable::index
(
    const List<scalar>& indices,
    const bool lastDim
) const
{
    label totalIndex = 0;

    forAll(dim_, i)
    {
        label dim = 1;
        for (int j = i + 1; j < dim_.size(); j++)
        {
            dim *= dim_[j] + 1;
        }

        totalIndex +=
            dim
           *Foam::min
            (
                Foam::max(label((indices[i] - min_[i])/delta_[i]), 0),
                dim_[i]
            );
    }

    if (lastDim)
    {
        label iLastdim = dim_.size() - 1;
        totalIndex += Foam::min
        (
            Foam::max
            (
                label((indices[iLastdim] - min_[iLastdim])/delta_[iLastdim]),
                0
            ),
            dim_[iLastdim]
        );
    }

    return totalIndex;
}


Foam::label
Foam::radiationModels::absorptionEmissionModels::interpolationLookUpTable::index
(
    const scalar indice
) const
{
    label i = 0;
    label totalIndex =
        Foam::min
        (
            Foam::max
            (
                label((indice - min_[i])/delta_[i]),
                0
            ),
            dim_[i]
        );

    return totalIndex;
}


bool
Foam::radiationModels::absorptionEmissionModels::interpolationLookUpTable::
checkRange
(
    const scalar lookUpValue,
    const label interfield
) const
{
    return lookUpValue >= min_[interfield] && lookUpValue <= max_[interfield];
}


Foam::scalar
Foam::radiationModels::absorptionEmissionModels::interpolationLookUpTable::
interpolate
(
    const label lo,
    const label hi,
    const scalar lookUpValue,
    const label ofield,
    const label interfield
) const
{
        if
        (
            List<scalarField>::operator[](interfield).operator[](hi)
         != List<scalarField>::operator[](interfield).operator[](lo)
        )
        {
            scalar output
            (
                List<scalarField>::operator[](ofield).operator[](lo)
              + (
                    List<scalarField>::operator[](ofield).operator[](hi)
                  - List<scalarField>::operator[](ofield).operator[](lo)
                )
               *(
                    lookUpValue
                  - List<scalarField>::operator[](interfield).operator[](lo)
                )
               /(
                    List<scalarField>::operator[](interfield).operator[](hi)
                  - List<scalarField>::operator[](interfield).operator[](lo)
                )
            );
            return output;
        }
        else
        {
            return List<scalarField>::operator[](ofield).operator[](lo);
        }
}


void
Foam::radiationModels::absorptionEmissionModels::interpolationLookUpTable::
dimensionTable()
{
    dim_.setSize(entries_.size());
    min_.setSize(entries_.size());
    delta_.setSize(entries_.size());
    max_.setSize(entries_.size());
    entryIndices_.setSize(entries_.size());
    outputIndices_.setSize(output_.size());
    label index = 0;
    label tableDim = 1;

    forAll(entries_,i)
    {
        dim_[i] = entries_[i].template lookup<label>("N");
        max_[i] = entries_[i].template lookup<scalar>("max");
        min_[i] = entries_[i].template lookup<scalar>("min");
        delta_[i] = (max_[i] - min_[i])/dim_[i];
        tableDim *= dim_[i] + 1;
        fieldIndices_.insert(entries_[i].lookup("name"), index);
        entryIndices_[i] = index;
        index++;
    }

    forAll(output_,i)
    {
        fieldIndices_.insert(output_[i].lookup("name"), index);
        outputIndices_[i] = index;
        index++;
    }

    List<scalarField>& internal = *this;

    internal.setSize(entries_.size() + output_.size());

    interpOutput_.setSize(entries_.size() + output_.size());

    forAll(internal, i)
    {
        internal[i].setSize(tableDim);
    }
}


void
Foam::radiationModels::absorptionEmissionModels::interpolationLookUpTable::
readTable
(
    const word& instance,
    const objectRegistry& obr
)
{
    IOdictionary control
    (
        IOobject
        (
            fileName_,
            instance,
            obr,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    control.lookup("fields") >> entries_;
    control.lookup("output") >> output_;
    control.lookup("values") >> *this;

    dimensionTable();

    check();

    if (this->size() == 0)
    {
        FatalErrorInFunction
            << "table is empty" << nl << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiationModels::absorptionEmissionModels::interpolationLookUpTable::
interpolationLookUpTable()
:
    List<scalarField>(),
    fileName_("fileNameIsUndefined")
{}


Foam::radiationModels::absorptionEmissionModels::interpolationLookUpTable::
interpolationLookUpTable
(
    const fileName& fn,
    const word& instance,
    const objectRegistry& obr
)
:
    List<scalarField>(),
    fileName_(fn),
    dim_(0),
    min_(0),
    delta_(0.0),
    max_(0.0),
    entries_(0),
    output_(0),
    entryIndices_(0),
    outputIndices_(0),
    interpOutput_(0)
{
    readTable(instance, obr);
}


Foam::radiationModels::absorptionEmissionModels::interpolationLookUpTable::
interpolationLookUpTable
(
     const interpolationLookUpTable& interpTable
)
:
    List<scalarField>(interpTable),
    fileName_(interpTable.fileName_),
    dim_(interpTable.dim_),
    min_(interpTable.min_),
    delta_(interpTable.delta_),
    max_(interpTable.max_),
    entries_(0),
    output_(0),
    entryIndices_(interpTable.entryIndices_),
    outputIndices_(interpTable.outputIndices_),
    interpOutput_(interpTable.interpOutput_)
{}


Foam::radiationModels::absorptionEmissionModels::interpolationLookUpTable::
interpolationLookUpTable
(
    const dictionary& dict
)
:
    List<scalarField>(),
    fileName_(fileName(dict.lookup("file")).expand()),
    dim_(0),
    min_(0.0),
    delta_(0.0),
    max_(0.0),
    entries_(dict.lookup("fields")),
    output_(dict.lookup("output")),
    entryIndices_(0),
    outputIndices_(0),
    fieldIndices_(0),
    interpOutput_(0)
{
    dimensionTable();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void
Foam::radiationModels::absorptionEmissionModels::interpolationLookUpTable::
check() const
{
    // check order in the first dimension.
    scalar prevValue = List<scalarField>::operator[](0).operator[](0);
    label dim = 1;
    for (int j = 1; j < dim_.size(); j++)
    {
        dim *= dim_[j] + 1;
    }

    for (label i = 1; i < dim_[0]; i++)
    {
        label index = i*dim;
        const scalar currValue =
            List<scalarField>::operator[](0).operator[](index);

        // avoid duplicate values (divide-by-zero error)
        if (currValue <= prevValue)
        {
            FatalErrorInFunction
                << "out-of-order value: " << currValue
                << " at index " << index << nl << exit(FatalError);
        }
        prevValue = currValue;
    }
}


void
Foam::radiationModels::absorptionEmissionModels::interpolationLookUpTable::
write
(
    Ostream& os,
    const fileName& fn,
    const word& instance,
    const objectRegistry& obr
) const
{
    IOdictionary control
    (
        IOobject
        (
            fn,
            instance,
            obr,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    );

    control.writeHeader(os);

    writeEntry(os, "fields", entries_);

    writeEntry(os, "output", output_);

    if (this->size() == 0)
    {
        FatalErrorInFunction
            << "table is empty" << nl << exit(FatalError);
    }
    writeEntry(os, "values", *this);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::scalarField&
Foam::radiationModels::absorptionEmissionModels::interpolationLookUpTable::
operator[]
(
    const label i
)
{
    const label n = this->size();

    if (n <= 1)
    {
        FatalErrorInFunction
            << "table has (" << n << ") columns" << nl << exit(FatalError);
    }
    else if (i < 0)
    {
        FatalErrorInFunction
            << "index (" << i << ") underflow" << nl << exit(FatalError);
    }
    else if (i >= n)
    {
        FatalErrorInFunction
            << "index (" << i << ") overflow" << nl << exit(FatalError);
    }

    return List<scalarField>::operator[](i);
}


const Foam::scalarField&
Foam::radiationModels::absorptionEmissionModels::interpolationLookUpTable::
operator[]
(
    const label i
) const
{
    const label n = this->size();

    if (n <= 1)
    {
        FatalErrorInFunction
            << "table has (" << n << ") columns" << nl << exit(FatalError);
    }
    else if (i < 0)
    {
        FatalErrorInFunction
            << "index (" << i << ") underflow" << nl << exit(FatalError);
    }
    else if (i >= n)
    {
        FatalErrorInFunction
            << "index (" << i << ") overflow" << nl
            << exit(FatalError);
    }

    return List<scalarField>::operator[](i);
}


bool
Foam::radiationModels::absorptionEmissionModels::interpolationLookUpTable::
found
(
    const word& fieldName
) const
{
    return fieldIndices_.found(fieldName);
}


const Foam::scalarList&
Foam::radiationModels::absorptionEmissionModels::interpolationLookUpTable::
lookUp
(
    const scalar retvals
)
{
    const label lo = index(retvals);
    findHi(lo, retvals);
    return interpOutput_;
}


void
Foam::radiationModels::absorptionEmissionModels::interpolationLookUpTable::
findHi
(
    const label lo,
    const scalar retvals
)
{
    forAll(outputIndices_,j)
    {
        scalar tmp = 0;
        label ofield = outputIndices_[j];
        scalar baseValue = List<scalarField>::operator[](ofield).operator[](lo);

        forAll(entryIndices_,i)
        {
            if (checkRange(retvals, entryIndices_[i]))
            {
                label dim = 1;

                label hi = Foam::min(lo + dim, (*this)[0].size() - 1);

                tmp += interpolate(lo, hi, retvals, ofield, entryIndices_[i])
                     - baseValue;
            }
            interpOutput_[entryIndices_[i]] = retvals;
        }

        tmp += baseValue;
        interpOutput_[outputIndices_[j]] = tmp;
    }
}


// ************************************************************************* //
