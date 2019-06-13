/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2019 OpenFOAM Foundation
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
    Foam::directPointPatchFieldMapper

Description
    direct pointPatchFieldMapper

\*---------------------------------------------------------------------------*/

#ifndef directPointPatchFieldMapper_H
#define directPointPatchFieldMapper_H

#include "pointPatchFieldMapper.H"
#include "directFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class directPointPatchFieldMapper Declaration
\*---------------------------------------------------------------------------*/

class directPointPatchFieldMapper
:
    public pointPatchFieldMapper,
    public directFieldMapper
{
    // Private Data

        //- Addressing from new back to old
        const labelUList& addressing_;

        //- Does map contain any unmapped values
        bool hasUnmapped_;


public:

    // Constructors

        //- Construct given addressing
        directPointPatchFieldMapper(const labelUList& addressing)
        :
            addressing_(addressing),
            hasUnmapped_(false)
        {
            if (addressing_.size() && min(addressing_) < 0)
            {
                hasUnmapped_ = true;
            }
        }


    //- Destructor
    virtual ~directPointPatchFieldMapper()
    {}


    // Member Functions

        bool hasUnmapped() const
        {
            return hasUnmapped_;
        }

        const labelUList& addressing() const
        {
            return addressing_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
