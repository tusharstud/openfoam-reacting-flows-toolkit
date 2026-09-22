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

\*---------------------------------------------------------------------------*/

#include "GurevichHSR.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::GurevichHSR<CloudType>::GurevichHSR
(
    const dictionary& dict,
    CloudType& owner
)
:
    SurfaceReactionModel<CloudType>(dict, owner, typeName),
    Ar_(this->coeffDict().template lookupOrDefault<scalar>("Ar", 1.5e4)),
    Ea_(this->coeffDict().template lookupOrDefault<scalar>("Ea", 8.372e4)),
    h_HSR_(this->coeffDict().template lookupOrDefault<scalar>("h_HSR", 3.1e7)),
    TmeltAl_(this->coeffDict().template lookupOrDefault<scalar>("TmeltAl", 933.0)),
    rhoAl_(this->coeffDict().template lookupOrDefault<scalar>("rhoAl", 2700.0)),
    AlSolidId_(-1),
    Al2O3SolidId_(-1),
    AlLiquidId_(-1),
    Al2O3LiquidId_(-1),
    O2GlobalId_(-1)
{
    // 1. Carrier O2 ID
    O2GlobalId_ = owner.composition().carrierId("O2", true);

    // 2. Solid components (Al core and Al2O3 layer)
    const label idSolid = owner.composition().idSolid();
    if (idSolid >= 0)
    {
        AlSolidId_ = owner.composition().localId(idSolid, "Al", true);
        if (AlSolidId_ < 0) AlSolidId_ = owner.composition().localId(idSolid, "C", true);
        if (AlSolidId_ < 0) AlSolidId_ = owner.composition().localId(idSolid, "aluminum", true);

        Al2O3SolidId_ = owner.composition().localId(idSolid, "Al2O3", true);
        if (Al2O3SolidId_ < 0) Al2O3SolidId_ = owner.composition().localId(idSolid, "ash", true);
        if (Al2O3SolidId_ < 0) Al2O3SolidId_ = owner.composition().localId(idSolid, "alumina", true);

        const wordList& sNames = owner.composition().componentNames(idSolid);
        if (AlSolidId_ < 0 && sNames.size() > 0) AlSolidId_ = 0;
        if (Al2O3SolidId_ < 0 && sNames.size() > 1) Al2O3SolidId_ = 1;
    }

    // 3. Liquid metal and oxide components
    const label idLiquid = owner.composition().idLiquid();
    if (idLiquid >= 0)
    {
        AlLiquidId_ = owner.composition().localId(idLiquid, "Al_liquid", true);
        if (AlLiquidId_ < 0) AlLiquidId_ = owner.composition().localId(idLiquid, "Al", true);
        if (AlLiquidId_ < 0) AlLiquidId_ = owner.composition().localId(idLiquid, "aluminum", true);

        Al2O3LiquidId_ = owner.composition().localId(idLiquid, "Al2O3_liquid", true);
        if (Al2O3LiquidId_ < 0) Al2O3LiquidId_ = owner.composition().localId(idLiquid, "Al2O3", true);
        if (Al2O3LiquidId_ < 0) Al2O3LiquidId_ = owner.composition().localId(idLiquid, "ash", true);

        const wordList& lNames = owner.composition().componentNames(idLiquid);
        if (AlLiquidId_ < 0 && lNames.size() > 0) AlLiquidId_ = 0;
        if (Al2O3LiquidId_ < 0 && lNames.size() > 1) Al2O3LiquidId_ = 1;
    }

    Info<< "    GurevichHSR: Al(s) ID = " << AlSolidId_
        << ", Al2O3(s) ID = " << Al2O3SolidId_
        << ", Al(l) ID = " << AlLiquidId_
        << ", Al2O3(l) ID = " << Al2O3LiquidId_
        << ", O2(carrier) ID = " << O2GlobalId_ << endl;
}


template<class CloudType>
Foam::GurevichHSR<CloudType>::GurevichHSR
(
    const GurevichHSR<CloudType>& srm
)
:
    SurfaceReactionModel<CloudType>(srm),
    Ar_(srm.Ar_),
    Ea_(srm.Ea_),
    h_HSR_(srm.h_HSR_),
    TmeltAl_(srm.TmeltAl_),
    rhoAl_(srm.rhoAl_),
    AlSolidId_(srm.AlSolidId_),
    Al2O3SolidId_(srm.Al2O3SolidId_),
    AlLiquidId_(srm.AlLiquidId_),
    Al2O3LiquidId_(srm.Al2O3LiquidId_),
    O2GlobalId_(srm.O2GlobalId_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::GurevichHSR<CloudType>::~GurevichHSR()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::GurevichHSR<CloudType>::calculate
(
    const scalar dt,
    const label celli,
    const scalar d,
    const scalar T,
    const scalar Tc,
    const scalar pc,
    const scalar rhoc,
    const scalar mass,
    const scalarField& YGas,
    const scalarField& YLiquid,
    const scalarField& YSolid,
    const scalarField& YMixture,
    const scalar N,
    scalarField& dMassGas,
    scalarField& dMassLiquid,
    scalarField& dMassSolid,
    scalarField& dMassSRCarrier
) const
{
    // Thermal ignition barrier: HSR activates once Al melts (T >= 933 K or Y_Al(l) > 0)
    const bool isIgnited = (T >= TmeltAl_) || (AlLiquidId_ >= 0 && AlLiquidId_ < YLiquid.size() && YLiquid[AlLiquidId_] > 1e-3);
    if (!isIgnited)
    {
        return 0.0;
    }

    const label idSolid = this->owner().composition().idSolid();
    const label idLiquid = this->owner().composition().idLiquid();

    scalar mAlSolid = 0.0;
    if (idSolid >= 0 && AlSolidId_ >= 0 && AlSolidId_ < YSolid.size())
    {
        mAlSolid = mass * YMixture[idSolid] * YSolid[AlSolidId_];
    }

    scalar mAlLiquid = 0.0;
    if (idLiquid >= 0 && AlLiquidId_ >= 0 && AlLiquidId_ < YLiquid.size())
    {
        mAlLiquid = mass * YMixture[idLiquid] * YLiquid[AlLiquidId_];
    }

    const scalar mAlTotal = mAlSolid + mAlLiquid;
    if (mAlTotal <= 1e-25)
    {
        return 0.0; // Aluminum fully consumed
    }

    // Effective diameter of aluminum core: dp,Al = (6*mp,Al / (pi*rhoAl))^(1/3) (Eq. 12)
    const scalar dpAl = Foam::cbrt(6.0 * mAlTotal / (Foam::constant::mathematical::pi * rhoAl_));
    const scalar ApEff = Foam::constant::mathematical::pi * sqr(dpAl);

    // Fuel vapor displacement of surface oxygen (Paper lines 560-565):
    // As Tp rises towards boiling point, Al vapor pressure blanket displaces surface oxidizer
    const scalar exponent = 12.19 - 34037.0 / max(T, scalar(300.0));
    const scalar pSat = 101325.0 * Foam::exp(min(exponent, scalar(15.0)));
    const scalar XFs = min(pSat / max(pc, scalar(1e4)), scalar(0.999));

    // Oxygen mass fraction at particle surface free of vapor blockage
    scalar Yox_s = 0.233 * (1.0 - XFs);
    if (O2GlobalId_ >= 0)
    {
        Yox_s = this->owner().composition().carrier().Y()[O2GlobalId_][celli] * (1.0 - XFs);
    }
    Yox_s = max(Yox_s, scalar(0.0));

    if (Yox_s <= 1e-6)
    {
        return 0.0;
    }

    // Universal gas constant: R_gas = 8.314462 J/(mol K)
    // Ea_ = 83720.0 J/mol (83.72 kJ/mol)
    const scalar R_gas = 8.314462;
    const scalar expTerm = Foam::exp(-Ea_ / (R_gas * max(T, scalar(300.0))));

    // Al consumption rate by HSR [kg/s] (Eq. 10)
    const scalar mDotHSR_Al = ApEff * rhoc * Yox_s * Ar_ * expTerm;
    scalar dMassAl = mDotHSR_Al * dt;
    dMassAl = min(dMassAl, mAlTotal);

    if (dMassAl <= 0.0)
    {
        return 0.0;
    }

    // Stoichiometry: Al + 0.75 O2 -> 0.5 Al2O3
    // W_Al = 26.98, W_O2 = 32.00, W_Al2O3 = 101.96
    const scalar stoichO2 = 0.75 * 32.00 / 26.98;     // ~ 0.88955 kg O2 / kg Al
    const scalar stoichAl2O3 = 0.5 * 101.96 / 26.98; // ~ 1.88955 kg Al2O3 / kg Al

    const scalar dMassO2 = stoichO2 * dMassAl;
    const scalar dMassAl2O3 = stoichAl2O3 * dMassAl;

    // In OpenFOAM convention:
    // massNew = mass0 - sum(dMass)
    // Hence, consumed mass is POSITIVE (+), and deposited/added mass is NEGATIVE (-)

    // 1. Consume Al from particle (liquid preferred if available, then solid)
    if (mAlLiquid > 0.0 && AlLiquidId_ >= 0 && AlLiquidId_ < dMassLiquid.size())
    {
        const scalar dLiq = min(dMassAl, mAlLiquid);
        dMassLiquid[AlLiquidId_] += dLiq;
        const scalar rem = dMassAl - dLiq;
        if (rem > 0.0 && AlSolidId_ >= 0 && AlSolidId_ < dMassSolid.size())
        {
            dMassSolid[AlSolidId_] += rem;
        }
    }
    else if (AlSolidId_ >= 0 && AlSolidId_ < dMassSolid.size())
    {
        dMassSolid[AlSolidId_] += dMassAl;
    }

    // 2. Deposit condensed Al2O3 onto particle
    // If oxide has already melted (T >= 2327 K) and liquid oxide exists, deposit to liquid oxide;
    // Otherwise deposit into solid oxide
    if (T >= 2327.0 && Al2O3LiquidId_ >= 0 && Al2O3LiquidId_ < dMassLiquid.size())
    {
        dMassLiquid[Al2O3LiquidId_] -= dMassAl2O3;
    }
    else if (Al2O3SolidId_ >= 0 && Al2O3SolidId_ < dMassSolid.size())
    {
        dMassSolid[Al2O3SolidId_] -= dMassAl2O3;
    }

    // 3. Deplete O2 from Eulerian carrier phase (negative source in carrier)
    if (O2GlobalId_ >= 0 && O2GlobalId_ < dMassSRCarrier.size())
    {
        dMassSRCarrier[O2GlobalId_] -= dMassO2;
    }

    // 4. Reaction heat release [J] (Eq. 14)
    const scalar hReaction = dMassAl * h_HSR_;

    return hReaction;
}

// ************************************************************************* //
