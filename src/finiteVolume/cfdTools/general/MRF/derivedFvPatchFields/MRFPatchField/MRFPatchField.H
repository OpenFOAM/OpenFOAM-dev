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

Class
    Foam::MRFnoSlipFvPatchVectorField

Description

See also
    Foam::MRFZone

SourceFiles
    MRFPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef MRFPatchField_H
#define MRFPatchField_H

#include "fvPatchFields.H"
#include "MRFZone.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class MRFPatchField Declaration
\*---------------------------------------------------------------------------*/

class MRFPatchField
{
    // Private data

        word MRFZoneName_;


public:

    // Constructors

        //- Construct null
        MRFPatchField();

        //- Construct from dictionary
        MRFPatchField(const dictionary&);


    // Member Functions

        const MRFZone& MRFzone(const objectRegistry& obr) const;

        void makeAbsolute(fvPatchField<vector>& Up) const;

        void makeRelative(fvPatchField<vector>& Up) const;

        void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
