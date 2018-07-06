/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "ensightPart.H"
#include "dictionary.H"
#include "ListOps.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
    defineTypeNameAndDebug(ensightPart, 0);
    defineTemplateTypeNameAndDebug(IOPtrList<ensightPart>, 0);
    defineRunTimeSelectionTable(ensightPart, istream);
}

const Foam::List<Foam::word> Foam::ensightPart::elemTypes_(0);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::ensightPart::isFieldDefined(const List<scalar>& field) const
{
    forAll(elemLists_, elemI)
    {
        const labelUList& idList = elemLists_[elemI];

        forAll(idList, i)
        {
            const label id = idList[i];

            if (id >= field.size() || std::isnan(field[id]))
            {
                return false;
            }
        }
    }
    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ensightPart::ensightPart
()
:
    number_(0),
    name_(""),
    elemLists_(0),
    offset_(0),
    size_(0),
    isCellData_(true),
    matId_(0),
    points_(pointField::null())
{}


Foam::ensightPart::ensightPart
(
    label partNumber,
    const string& partDescription
)
:
    number_(partNumber),
    name_(partDescription),
    elemLists_(0),
    offset_(0),
    size_(0),
    isCellData_(true),
    matId_(0),
    points_(pointField::null())
{}


Foam::ensightPart::ensightPart
(
    label partNumber,
    const string& partDescription,
    const pointField& points
)
:
    number_(partNumber),
    name_(partDescription),
    elemLists_(0),
    offset_(0),
    size_(0),
    isCellData_(true),
    matId_(0),
    points_(points)
{}


Foam::ensightPart::ensightPart(const ensightPart& part)
:
    number_(part.number_),
    name_(part.name_),
    elemLists_(part.elemLists_),
    offset_(part.offset_),
    size_(part.size_),
    isCellData_(part.isCellData_),
    matId_(part.matId_),
    points_(part.points_)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::ensightPart> Foam::ensightPart::New(Istream& is)
{
    const word partType(is);

    istreamConstructorTable::iterator cstrIter =
        istreamConstructorTablePtr_->find(partType);

    if (cstrIter == istreamConstructorTablePtr_->end())
    {
        FatalIOErrorInFunction
        (
            is
        )   << "unknown ensightPart type "
            << partType << nl << nl
            << "Valid ensightPart types are :" << endl
            << istreamConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return autoPtr<ensightPart>(cstrIter()(is));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ensightPart::~ensightPart()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ensightPart::renumber(const labelUList& origId)
{
    // transform to global values first
    if (offset_)
    {
        forAll(elemLists_, elemI)
        {
            labelList& idList = elemLists_[elemI];
            forAll(idList, i)
            {
                idList[i] += offset_;
            }
        }

        offset_ = 0;
    }

    if (origId.size())
    {
        forAll(elemLists_, elemI)
        {
            inplaceRenumber(origId, elemLists_[elemI]);
        }
    }
}


// ************************************************************************* //
