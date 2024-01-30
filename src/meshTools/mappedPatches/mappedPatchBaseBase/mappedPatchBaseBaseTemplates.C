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

#include "mappedPatchBaseBase.H"
#include "stringOps.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class PatchFieldType, class FieldType>
void Foam::mappedPatchBaseBase::validateMapForField
(
    const PatchFieldType& field,
    const FieldType& iF,
    const dictionary& context,
    const label froms
)
{
    const polyPatch& pp = field.patch().patch();

    if (!isA<mappedPatchBaseBase>(pp))
    {
        OStringStream str;
        str << "Field " << iF.name() << " of type "
            << field.type() << " cannot apply to patch " << pp.name()
            << " because the patch is not of " << typeName << " type";
        FatalIOErrorInFunction(context)
            << stringOps::breakIntoIndentedLines(str.str()).c_str()
            << exit(FatalIOError);
    }

    refCast<const mappedPatchBaseBase>(pp).validateForField
    (
        field,
        iF,
        context,
        froms
    );
}


template<class PatchFieldType, class FieldType>
void Foam::mappedPatchBaseBase::validateForField
(
    const PatchFieldType& field,
    const FieldType& iF,
    const dictionary& context,
    const label froms
) const
{
    const bool isNotRegion = !sameRegion() && (froms & from::sameRegion);
    const bool isRegion = sameRegion() && (froms & from::differentRegion);
    const bool isPatch = samePatch() && (froms & from::differentPatch);

    OStringStream str;

    if (isNotRegion || isRegion || isPatch)
    {
        str << "Field " << iF.name() << " of type "
            << field.type() << " cannot apply to patch " << patch_.name()
            << " because values are mapped from ";
    }

    if (isNotRegion)
    {
        str << "a different region";
    }
    else if (isRegion)
    {
        str << "within the same region";
    }
    else if (isPatch)
    {
        str << "the same patch";
    }

    if (isNotRegion || isRegion || isPatch)
    {
        FatalIOErrorInFunction(context)
            << stringOps::breakIntoIndentedLines(str.str()).c_str()
            << exit(FatalIOError);
    }
}


// ************************************************************************* //
