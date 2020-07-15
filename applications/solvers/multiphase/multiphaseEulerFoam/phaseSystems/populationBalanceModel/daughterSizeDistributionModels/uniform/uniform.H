/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2020 OpenFOAM Foundation
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
    Foam::diameterModels::daughterSizeDistributionModels::uniform

Description
    Uniform daughter size distribution.

SourceFiles
    uniform.C

\*---------------------------------------------------------------------------*/

#ifndef uniform_H
#define uniform_H

#include "daughterSizeDistributionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace daughterSizeDistributionModels
{

/*---------------------------------------------------------------------------*\
                           Class uniform Declaration
\*---------------------------------------------------------------------------*/

class uniform
:
    public daughterSizeDistributionModel
{
public:

    //- Runtime type information
    TypeName("uniform");


    // Constructor

        uniform
        (
            const breakupModel& breakup,
            const dictionary& dict
        );


    //- Destructor
    virtual ~uniform();


    // Member Functions

        //- Calculate contribution to sizeGroup i due to breakup in sizeGroup k
        virtual dimensionedScalar calcNik(const label i, const label k) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace daughterSizeDistributionModels
} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
