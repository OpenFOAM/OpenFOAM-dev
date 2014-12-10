/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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

#include "CodedSource.H"
#include "stringOps.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::fv::CodedSource<Type>::writeData(Ostream& os) const
{
    os  << indent << name_ << endl;
    dict_.write(os);
}


template<class Type>
bool Foam::fv::CodedSource<Type>::read(const dictionary& dict)
{
    if (option::read(dict))
    {
        coeffs_.lookup("fieldNames") >> fieldNames_;
        applied_.setSize(fieldNames_.size(), false);
        coeffs_.lookup("redirectType") >> redirectType_;

        // Code snippets
        {
            const entry& e = coeffs_.lookupEntry
            (
                "codeCorrect",
                false,
                false
            );
            codeCorrect_ = stringOps::trim(e.stream());
            stringOps::inplaceExpand(codeCorrect_, coeffs_);
            dynamicCodeContext::addLineDirective
            (
                codeCorrect_,
                e.startLineNumber(),
                coeffs_.name()
            );
        }

        {
            const entry& e = coeffs_.lookupEntry
            (
                "codeAddSup",
                false,
                false
            );
            codeAddSup_ = stringOps::trim(e.stream());
            stringOps::inplaceExpand(codeAddSup_, coeffs_);
            dynamicCodeContext::addLineDirective
            (
                codeAddSup_,
                e.startLineNumber(),
                coeffs_.name()
            );
        }

        {
            const entry& e = coeffs_.lookupEntry
            (
                "codeSetValue",
                false,
                false
            );
            codeSetValue_ = stringOps::trim(e.stream());
            stringOps::inplaceExpand(codeSetValue_, coeffs_);
            dynamicCodeContext::addLineDirective
            (
                codeSetValue_,
                e.startLineNumber(),
                coeffs_.name()
            );
        }

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
