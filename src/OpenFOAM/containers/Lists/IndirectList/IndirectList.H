/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::IndirectList

Description
    A List with indirect addressing.

See also
    Foam::UIndirectList for a version without any allocation for the
    addressing.

SourceFiles
    IndirectListI.H

\*---------------------------------------------------------------------------*/

#ifndef IndirectList_H
#define IndirectList_H

#include "UIndirectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class IndirectListAddressing Declaration
\*---------------------------------------------------------------------------*/

//- A helper class for storing addresses.
class IndirectListAddressing
{
    // Private data

        //- Storage for the list addressing
        List<label> addressing_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        IndirectListAddressing(const IndirectListAddressing&);

        //- Disallow default bitwise assignment
        void operator=(const IndirectListAddressing&);


protected:

    // Constructors

        //- Construct by copying the addressing array
        explicit inline IndirectListAddressing(const labelUList& addr);

        //- Construct by transferring addressing array
        explicit inline IndirectListAddressing(const Xfer<List<label>>& addr);


    // Member Functions

        // Access

            //- Return the list addressing
            inline const List<label>& addressing() const;

        // Edit

            //- Reset addressing
            inline void resetAddressing(const labelUList&);
            inline void resetAddressing(const Xfer<List<label>>&);

};


/*---------------------------------------------------------------------------*\
                        Class IndirectList Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class IndirectList
:
    private IndirectListAddressing,
    public  UIndirectList<T>
{
    // Private Member Functions

        //- Disallow default assignment operator
        void operator=(const IndirectList<T>&);

        //- Disallow assignment from UIndirectList
        void operator=(const UIndirectList<T>&);


public:

    // Constructors

        //- Construct given the complete list and the addressing array
        inline IndirectList(const UList<T>&, const labelUList&);

        //- Construct given the complete list and by transferring addressing
        inline IndirectList(const UList<T>&, const Xfer<List<label>>&);

        //- Copy constructor
        inline IndirectList(const IndirectList<T>&);

        //- Construct from UIndirectList
        explicit inline IndirectList(const UIndirectList<T>&);


    // Member Functions


        // Access

            //- Return the list addressing
            using UIndirectList<T>::addressing;


        // Edit

            //- Reset addressing
            using IndirectListAddressing::resetAddressing;


        // Member Operators

            //- Assignment operator
            using UIndirectList<T>::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "IndirectListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
