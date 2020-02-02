/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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

#include "CodedBase.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class CodedType>
const Foam::word Foam::CodedBase<CodedType>::codeTemplateC =
    CodedType::typeName + "Template.C";

template<class CodedType>
const Foam::word Foam::CodedBase<CodedType>::codeTemplateH =
    CodedType::typeName + "Template.H";


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CodedType>
const Foam::wordList&
Foam::CodedBase<CodedType>::codeKeys() const
{
    return codeKeys_;
}


template<class CodedType>
Foam::string Foam::CodedBase<CodedType>::description() const
{
    return CodedType::typeName + " " + codeName();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CodedType>
Foam::CodedBase<CodedType>::CodedBase(const dictionary& dict)
:
    codeName_(dict.lookup("name")),
    dict_(dict)
{}


template<class CodedType>
Foam::CodedBase<CodedType>::CodedBase(const CodedBase<CodedType>& cb)
:
    codeName_(cb.codeName_),
    dict_(cb.dict_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


template<class CodedType>
Foam::CodedBase<CodedType>::~CodedBase()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CodedType>
const Foam::word& Foam::CodedBase<CodedType>::codeName() const
{
    return codeName_;
}


template<class CodedType>
const Foam::dictionary&
Foam::CodedBase<CodedType>::codeDict() const
{
    return dict_;
}


template<class CodedType>
void Foam::CodedBase<CodedType>::writeCode(Ostream& os) const
{
    writeEntry(os, "name", codeName_);

    forAll(codeKeys_, i)
    {
        if (dict_.found(codeKeys_[i]))
        {
            writeKeyword(os, codeKeys_[i]);
            os.write(verbatimString(dict_[codeKeys_[i]]))
            << token::END_STATEMENT << nl;
        }
    }
}


// ************************************************************************* //
