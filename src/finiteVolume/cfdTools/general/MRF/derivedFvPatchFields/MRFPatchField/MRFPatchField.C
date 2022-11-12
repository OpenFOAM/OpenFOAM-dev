/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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

#include "MRFPatchField.H"
#include "IOMRFZoneList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MRFPatchField::MRFPatchField()
{}


Foam::MRFPatchField::MRFPatchField(const dictionary& dict)
:
    MRFZoneName_(dict.lookupOrDefault("MRFZoneName", word::null))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::MRFZone& Foam::MRFPatchField::MRFzone
(
    const objectRegistry& obr
) const
{
    // Get reference to the MRF model
    const MRFZoneList& mrf =
        obr.lookupObject<IOMRFZoneList>("MRFProperties");

    if (MRFZoneName_ != word::null)
    {
        forAll(mrf, i)
        {
            if (mrf[i].name() == MRFZoneName_)
            {
                return mrf[i];
            }
        }

        FatalErrorInFunction
            << "Cannot find MRF zone " << MRFZoneName_
            << exit(FatalError);
    }
    else if (mrf.size() == 1)
    {
        return mrf[0];
    }
    else if (mrf.size() == 0)
    {
        FatalErrorInFunction
            << "There are no MRF zones"
            << exit(FatalError);
    }
    else
    {
        FatalErrorInFunction
            << "MRFZoneName not specified"
            << exit(FatalError);
    }

    return mrf[0];
}


void Foam::MRFPatchField::makeAbsolute(fvPatchField<vector>& Up) const
{
    MRFzone(Up.db()).makeAbsolute(Up, Up.patch().index());
}


void Foam::MRFPatchField::makeRelative(fvPatchField<vector>& Up) const
{
    MRFzone(Up.db()).makeRelative(Up, Up.patch().index());
}


void Foam::MRFPatchField::write(Ostream& os) const
{
    if (MRFZoneName_ != word::null)
    {
        writeEntry(os, "MRFZoneName", MRFZoneName_);
    }
}


// ************************************************************************* //
