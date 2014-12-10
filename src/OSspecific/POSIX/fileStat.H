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
    Foam::fileStat

Description
    Wrapper for stat() system call.

Warning
    on Linux (an maybe on others) a stat() of an nfs mounted (remote)
    file does never timeout and cannot be interrupted!
    So e.g. Foam::ping first and hope nfs is running.

SourceFiles
    fileStat.C

\*---------------------------------------------------------------------------*/

#ifndef fileStat_H
#define fileStat_H

#include <sys/stat.h>
#include <sys/types.h>

#include "label.H"
#include "fileName.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class fileStat;

Istream& operator>>(Istream&, fileStat&);
Ostream& operator<<(Ostream&, const fileStat&);


/*---------------------------------------------------------------------------*\
                           Class fileStat Declaration
\*---------------------------------------------------------------------------*/

class fileStat
{
    // Private data

        struct stat status_;

        bool isValid_;


public:

    // Constructors

        //- Empty constructor
        fileStat();

        //- Construct from components
        fileStat(const fileName& fName, const unsigned int maxTime=0);

        //- Construct from Istream
        fileStat(Istream&);


    // Member Functions

        // Access

            //- Raw status
            const struct stat& status() const
            {
                return status_;
            }

            //- Did constructor fail
            bool isValid() const
            {
                return isValid_;
            }


        // Check

            //- compare two fileStats for same device
            bool sameDevice(const fileStat& stat2) const;

            //- compare two fileStats for same Inode
            bool sameINode(const fileStat& stat2) const;

            //- compare state against inode
            bool sameINode(const label iNode) const;


    // IOstream Operators

        friend Istream& operator>>(Istream&, fileStat&);
        friend Ostream& operator<<(Ostream&, const fileStat&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
