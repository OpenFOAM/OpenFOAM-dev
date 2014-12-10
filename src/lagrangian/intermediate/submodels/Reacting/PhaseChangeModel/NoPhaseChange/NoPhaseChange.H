/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::NoPhaseChange

Description
    Dummy phase change model for 'none'

\*---------------------------------------------------------------------------*/

#ifndef NoPhaseChange_H
#define NoPhaseChange_H

#include "PhaseChangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                       Class NoPhaseChange Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NoPhaseChange
:
    public PhaseChangeModel<CloudType>
{
public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from dictionary
        NoPhaseChange(const dictionary&, CloudType&);

        //- Construct copy
        NoPhaseChange(const NoPhaseChange<CloudType>& pcm);

        //- Construct and return a clone
        virtual autoPtr<PhaseChangeModel<CloudType> > clone() const
        {
            return autoPtr<PhaseChangeModel<CloudType> >
            (
                new NoPhaseChange<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~NoPhaseChange();


    // Member Functions

        //- Flag to indicate whether model activates phase change model
        virtual bool active() const;

        //- Update model
        virtual void calculate
        (
            const scalar dt,
            const label cellI,
            const scalar Re,
            const scalar d,
            const scalar nu,
            const scalar T,
            const scalar Ts,
            const scalar pc,
            scalarField& dMassPC
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "NoPhaseChange.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
