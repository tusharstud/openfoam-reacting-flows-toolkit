/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::Kavanau

Description
    Kavanau correlation for transition-regime heat transfer of micron-sized
    particles (Zhang et al., 2022, Eq. 8).

\*---------------------------------------------------------------------------*/

#ifndef Kavanau_H
#define Kavanau_H

#include "HeatTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class Kavanau Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class Kavanau
:
    public HeatTransferModel<CloudType>
{

public:

    //- Runtime type information
    TypeName("Kavanau");


    // Constructors

        //- Construct from dictionary
        Kavanau(const dictionary& dict, CloudType& cloud);

        //- Construct copy
        Kavanau(const Kavanau<CloudType>& im);

        //- Construct and return a clone
        virtual autoPtr<HeatTransferModel<CloudType>> clone() const
        {
            return autoPtr<HeatTransferModel<CloudType>>
            (
                new Kavanau<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~Kavanau();


    // Member Functions

        // Evaluation

            //- Nusselt number
            virtual scalar Nu
            (
                const scalar Re,
                const scalar Pr
            ) const;

            //- Heat transfer coefficient
            virtual scalar htc
            (
                const scalar dp,
                const scalar Re,
                const scalar Pr,
                const scalar kappa,
                const scalar NCpW
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Kavanau.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
