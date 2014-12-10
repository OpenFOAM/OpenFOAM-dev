/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::TimeScaleModels::nonEquilibrium

Description
    Non-Equlibrium model for the time scale over which properties of a
    dispersed phase tend towards the mean value.

    Reference:
    \verbatim
        "An improved collision damping time for MP-PIC calculations of dense
        particle flows with applications to polydisperse sedimenting beds and
        colliding particle jets"
        P O'Rourke and D Snider
        Chemical Engineering Science
        Volume 65, Issue 22, Pages 6014-6028, November 2010
    \endverbatim

SourceFiles
    nonEquilibrium.C

\*---------------------------------------------------------------------------*/

#ifndef nonEquilibrium_H
#define nonEquilibrium_H

#include "TimeScaleModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace TimeScaleModels
{

/*---------------------------------------------------------------------------*\
                        Class nonEquilibrium Declaration
\*---------------------------------------------------------------------------*/

class nonEquilibrium
:
    public TimeScaleModel
{
public:

    //- Runtime type information
    TypeName("nonEquilibrium");


    //- Constructors

        //- Construct from components
        nonEquilibrium(const dictionary& dict);

        //- Construct a copy
        nonEquilibrium(const nonEquilibrium& hc);

        //- Construct and return a clone
        virtual autoPtr<TimeScaleModel> clone() const
        {
            return autoPtr<TimeScaleModel>
            (
                new nonEquilibrium(*this)
            );
        }


    //- Destructor
    virtual ~nonEquilibrium();


    //- Member Functions

        //- Time scale
        tmp<FieldField<Field, scalar> > oneByTau
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
