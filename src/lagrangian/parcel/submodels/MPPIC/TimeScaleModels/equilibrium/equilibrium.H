/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2021 OpenFOAM Foundation
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
    Foam::TimeScaleModels::equilibrium

Description
    Equilibrium model for the time scale over which properties of a dispersed
    phase tend towards the mean value.

    Reference:
    \verbatim
        "A model for collisional exchange in gas/liquid/solid fluidised beds"
        P O'Rourke, P Zhao and D Snider
        Chemical Engineering Science
        Volume 64, Issue 8, Pages 1784-1797, April 2009
    \endverbatim

SourceFiles
    equilibrium.C

\*---------------------------------------------------------------------------*/

#ifndef equilibrium_H
#define equilibrium_H

#include "TimeScaleModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace TimeScaleModels
{

/*---------------------------------------------------------------------------*\
                         Class equilibrium Declaration
\*---------------------------------------------------------------------------*/

class equilibrium
:
    public TimeScaleModel
{
public:

    //- Runtime type information
    TypeName("equilibrium");


    //- Constructors

        //- Construct from components
        equilibrium(const dictionary& dict);

        //- Construct a copy
        equilibrium(const equilibrium& hc);

        //- Construct and return a clone
        virtual autoPtr<TimeScaleModel> clone() const
        {
            return autoPtr<TimeScaleModel>
            (
                new equilibrium(*this)
            );
        }


    //- Destructor
    virtual ~equilibrium();


    //- Member Functions

        //- Time scale
        tmp<FieldField<Field, scalar>> oneByTau
        (
            const FieldField<Field, scalar>& alpha,
            const FieldField<Field, scalar>& r32,
            const FieldField<Field, scalar>& uSqr,
            const FieldField<Field, scalar>& f
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace TimeScaleModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
