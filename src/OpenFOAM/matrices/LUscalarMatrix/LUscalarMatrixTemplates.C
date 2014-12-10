/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "LUscalarMatrix.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::LUscalarMatrix::solve(Field<Type>& sourceSol) const
{
    if (Pstream::parRun())
    {
        Field<Type> completeSourceSol(n());

        if (Pstream::master(comm_))
        {
            typename Field<Type>::subField
            (
                completeSourceSol,
                sourceSol.size()
            ).assign(sourceSol);

            for
            (
                int slave=Pstream::firstSlave();
                slave<=Pstream::lastSlave(comm_);
                slave++
            )
            {
                IPstream::read
                (
                    Pstream::scheduled,
                    slave,
                    reinterpret_cast<char*>
                    (
                        &(completeSourceSol[procOffsets_[slave]])
                    ),
                    (procOffsets_[slave+1]-procOffsets_[slave])*sizeof(Type),
                    Pstream::msgType(),
                    comm_
                );
            }
        }
        else
        {
            OPstream::write
            (
                Pstream::scheduled,
                Pstream::masterNo(),
                reinterpret_cast<const char*>(sourceSol.begin()),
                sourceSol.byteSize(),
                Pstream::msgType(),
                comm_
            );
        }

        if (Pstream::master(comm_))
        {
            LUBacksubstitute(*this, pivotIndices_, completeSourceSol);

            sourceSol = typename Field<Type>::subField
            (
                completeSourceSol,
                sourceSol.size()
            );

            for
            (
                int slave=Pstream::firstSlave();
                slave<=Pstream::lastSlave(comm_);
                slave++
            )
            {
                OPstream::write
                (
                    Pstream::scheduled,
                    slave,
                    reinterpret_cast<const char*>
                    (
                        &(completeSourceSol[procOffsets_[slave]])
                    ),
                    (procOffsets_[slave + 1]-procOffsets_[slave])*sizeof(Type),
                    Pstream::msgType(),
                    comm_
                );
            }
        }
        else
        {
            IPstream::read
            (
                Pstream::scheduled,
                Pstream::masterNo(),
                reinterpret_cast<char*>(sourceSol.begin()),
                sourceSol.byteSize(),
                Pstream::msgType(),
                comm_
            );
        }
    }
    else
    {
        LUBacksubstitute(*this, pivotIndices_, sourceSol);
    }
}


// ************************************************************************* //
