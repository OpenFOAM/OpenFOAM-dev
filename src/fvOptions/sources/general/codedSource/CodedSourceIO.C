/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2020 OpenFOAM Foundation
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
bool Foam::fv::CodedSource<Type>::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        coeffs_.lookup("fields") >> fieldNames_;
        applied_.setSize(fieldNames_.size(), false);

        // The name keyword is "name". "redirectType" is also maintained here
        // for backwards compatibility, but "name" is taken in preference and
        // is printed in the error message if neither keyword is present.
        name_ = word::null;
        name_ = dict.lookupOrDefault("redirectType", name_);
        name_ = dict.lookupOrDefault("name", name_);
        if (name_ == word::null)
        {
            dict.lookup("name"); // <-- generate error message with "name" in it
        }

        updateLibrary();

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
