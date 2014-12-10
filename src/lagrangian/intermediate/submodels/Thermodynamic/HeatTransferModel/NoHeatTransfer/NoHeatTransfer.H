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
    Foam::NoHeatTransfer

Description
    Dummy heat transfer model for 'none'

\*---------------------------------------------------------------------------*/

#ifndef NoHeatTransfer_H
#define NoHeatTransfer_H

#include "HeatTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class NoHeatTransfer Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NoHeatTransfer
:
    public HeatTransferModel<CloudType>
{
public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from dictionary
        NoHeatTransfer(const dictionary&, CloudType& owner);

        //- Construct copy
        NoHeatTransfer(const NoHeatTransfer<CloudType>& im);

        //- Construct and return a clone
        virtual autoPtr<HeatTransferModel<CloudType> > clone() const
        {
            return autoPtr<HeatTransferModel<CloudType> >
            (
                new NoHeatTransfer<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~NoHeatTransfer();


    // Member Functions

        //- Flag to indicate whether model activates heat transfer model
        virtual bool active() const;

        //- Nusselt number
        virtual scalar Nu(const scalar, const scalar) const;

        //- Prandtl number
        virtual scalar Pr() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "NoHeatTransfer.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
