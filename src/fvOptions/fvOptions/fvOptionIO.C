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

#include "fvOption.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::option::writeHeader(Ostream& os) const
{
    os  << indent << name_ << nl
        << indent << token::BEGIN_BLOCK << incrIndent << nl;
}


void Foam::fv::option::writeFooter(Ostream& os) const
{
    os  << decrIndent << indent << token::END_BLOCK << endl;
}


void Foam::fv::option::writeData(Ostream& os) const
{
    os.writeKeyword("active") << active_ << token::END_STATEMENT << nl;
    os.writeKeyword("timeStart") << timeStart_ << token::END_STATEMENT << nl;
    os.writeKeyword("duration") << duration_ << token::END_STATEMENT << nl;
    os.writeKeyword("selectionMode")
        << selectionModeTypeNames_[selectionMode_] << nl;

    switch (selectionMode_)
    {
        case smPoints:
        {
            os.writeKeyword("points") << points_
                << token::END_STATEMENT << nl;
            break;
        }
        case smCellSet:
        {
            os.writeKeyword("cellSet") << cellSetName_
                << token::END_STATEMENT << nl;
            break;
        }
        case smCellZone:
        {
            os.writeKeyword("cellZone") << cellSetName_
                << token::END_STATEMENT << nl;
            break;
        }
        case smAll:
        {
            break;
        }
        case smMapRegion:
        {
            break;
        }
        default:
        {
            FatalErrorIn("option::writeData(Ostream&) const")
                << "Unknown selectionMode "
                << selectionMode_
                << abort(FatalError);
        }
    }
}


bool Foam::fv::option::read(const dictionary& dict)
{
    active_ = readBool(dict.lookup("active"));

    if (dict.readIfPresent("timeStart", timeStart_))
    {
        dict.lookup("duration") >> duration_;
    }

    coeffs_ = dict.subDict(modelType_ + "Coeffs");

    return true;
}


// ************************************************************************* //
