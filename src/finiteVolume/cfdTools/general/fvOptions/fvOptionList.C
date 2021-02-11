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

#include "fvOptionList.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(optionList, 0);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

const Foam::dictionary& Foam::fv::optionList::optionsDict
(
    const dictionary& dict
) const
{
    if (dict.found("options"))
    {
        return dict.subDict("options");
    }
    else
    {
        return dict;
    }
}


bool Foam::fv::optionList::readOptions(const dictionary& dict)
{
    const dictionary& optionsDict(this->optionsDict(dict));

    checkTimeIndex_ = mesh_.time().timeIndex() + 1;

    bool allOk = true;
    forAll(*this, i)
    {
        option& bs = this->operator[](i);
        bool ok = bs.read(optionsDict.subDict(bs.name()));
        allOk = (allOk && ok);
    }
    return allOk;
}


void Foam::fv::optionList::checkApplied() const
{
    if (mesh_.time().timeIndex() > checkTimeIndex_)
    {
        forAll(*this, i)
        {
            const option& source = this->operator[](i);

            wordHashSet notAddSupFields(source.addSupFields());
            notAddSupFields -= addSupFields_[i];

            forAllConstIter(wordHashSet, notAddSupFields, iter)
            {
                WarningInFunction
                    << "Source " << source.name() << " defined for field "
                    << iter.key() << " but never used" << endl;
            }

            wordHashSet notConstrainedFields(source.constrainedFields());
            notConstrainedFields -= constrainedFields_[i];

            forAllConstIter(wordHashSet, notConstrainedFields, iter)
            {
                WarningInFunction
                    << "Constraint " << source.name() << " defined for field "
                    << iter.key() << " but never used" << endl;
            }

            wordHashSet notCorrectedFields(source.correctedFields());
            notCorrectedFields -= correctedFields_[i];

            forAllConstIter(wordHashSet, notCorrectedFields, iter)
            {
                WarningInFunction
                    << "Correction " << source.name() << " defined for field "
                    << iter.key() << " but never used" << endl;
            }
        }

        checkTimeIndex_ = mesh_.time().timeIndex();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::optionList::optionList(const fvMesh& mesh, const dictionary& dict)
:
    PtrListDictionary<option>(0),
    mesh_(mesh),
    checkTimeIndex_(mesh_.time().timeIndex() + 1),
    addSupFields_(),
    constrainedFields_(),
    correctedFields_()
{
    reset(dict);
}


Foam::fv::optionList::optionList(const fvMesh& mesh)
:
    PtrListDictionary<option>(0),
    mesh_(mesh),
    checkTimeIndex_(mesh_.time().timeIndex() + 1),
    addSupFields_(),
    constrainedFields_(),
    correctedFields_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::optionList::reset(const dictionary& dict)
{
    const dictionary& optionsDict(this->optionsDict(dict));

    // Count number of active fvOptions
    label count = 0;
    forAllConstIter(dictionary, optionsDict, iter)
    {
        if (iter().isDict())
        {
            count++;
        }
    }

    this->setSize(count);

    addSupFields_.setSize(count);
    constrainedFields_.setSize(count);
    correctedFields_.setSize(count);

    label i = 0;
    forAllConstIter(dictionary, optionsDict, iter)
    {
        if (iter().isDict())
        {
            const word& name = iter().keyword();
            const dictionary& sourceDict = iter().dict();

            this->set
            (
                i,
                name,
                option::New(name, sourceDict, mesh_).ptr()
            );

            addSupFields_.set(i, new wordHashSet());
            constrainedFields_.set(i, new wordHashSet());
            correctedFields_.set(i, new wordHashSet());

            i++;
        }
    }
}


bool Foam::fv::optionList::addsSupToField(const word& fieldName) const
{
    forAll(*this, i)
    {
        if (this->operator[](i).addsSupToField(fieldName))
        {
            return true;
        }
    }

    return false;
}


bool Foam::fv::optionList::constrainsField(const word& fieldName) const
{
    forAll(*this, i)
    {
        if (this->operator[](i).constrainsField(fieldName))
        {
            return true;
        }
    }

    return false;
}


bool Foam::fv::optionList::correctsField(const word& fieldName) const
{
    forAll(*this, i)
    {
        if (this->operator[](i).correctsField(fieldName))
        {
            return true;
        }
    }

    return false;
}


bool Foam::fv::optionList::read(const dictionary& dict)
{
    return readOptions(dict);
}


bool Foam::fv::optionList::writeData(Ostream& os) const
{
    // Write list contents
    forAll(*this, i)
    {
        os  << nl;
        this->operator[](i).writeHeader(os);
        this->operator[](i).writeData(os);
        this->operator[](i).writeFooter(os);
    }

    // Check state of IOstream
    return os.good();
}


void Foam::fv::optionList::updateMesh(const mapPolyMesh& mpm)
{
    forAll(*this, i)
    {
        this->operator[](i).updateMesh(mpm);
    }
}


bool Foam::fv::optionList::movePoints()
{
    bool allOk = true;

    forAll(*this, i)
    {
        allOk = allOk && this->operator[](i).movePoints();
    }

    return allOk;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const fv::optionList& options)
{
    options.writeData(os);
    return os;
}


// ************************************************************************* //
