!!****if* source/physics/Eos/EosMain/Eos_getData_loop1
!! NAME
!!
!!  Eos_getData_loop1
!! 
!! SYNOPSIS
!!
!!  call Eos_getData(  )
!!
!! DESCRIPTION
!! Inner loop that was in Eos_getData is outlined here
!!***

! solnData depends on the ordering on unk
!!REORDER(4): solnData

subroutine Eos_getData_loop1(vecLen, eosData, solnData, i,j,k,n, pres_map,dens_map,gamc_map,&
                    game_map,temp_map,entr_map,eint_map,ener_map, velx_map, vely_map, velz_map, &
                    sumy_map, ye_map, pres,dens,gamc,temp,abar,zbar,eint,ekin,entr)

use Eos_data, ONLY: eos_eintSwitch, eos_smalle

implicit none

#include "Eos.h"
#include "Eos_map.h"
#include "constants.h"
#include "Flash.h"

        integer, intent(in) :: vecLen
        real, dimension(EOS_NUM*vecLen),intent(INOUT) :: eosData
        real, pointer:: solnData(:,:,:,:)
        integer, intent(in) :: i,j,k
        integer, intent(inout) :: n
        integer, intent(in) :: pres,dens,gamc,temp,abar,zbar,eint,ekin,entr
        integer, intent(in) :: pres_map,dens_map,gamc_map,game_map,temp_map,entr_map
        integer, intent(in) :: eint_map,ener_map, velx_map, vely_map, velz_map, sumy_map, ye_map

        real :: kineticEnergy, internalEnergy
           if (velx_map > 0 .AND. vely_map > 0 .AND. velz_map > 0) then
              kineticEnergy  = 0.5*(solnData(velx_map,i,j,k)**2 + &
                                    solnData(vely_map,i,j,k)**2 + &
                                    solnData(velz_map,i,j,k)**2)
           else
              kineticEnergy = 0.0
           end if
           
           n=n+1
           eosData(ekin+n) = kineticEnergy
           !! kineticEnergy holds velocity vector information -- 1/2 * Vmag**2
           !! internalEnergy holds eint (directly)  or energyTotal - ekinetic (calculated),
           !!          depending upon eintSwitch
           if(eint_map /= NONEXISTENT) then
              internalEnergy  = solnData(eint_map,i,j,k)
              if(ener_map /= NONEXISTENT) then
                 if ( solnData(ener_map,i,j,k) - kineticEnergy > max(eos_smalle, eos_eintSwitch*kineticEnergy)) then
                    internalEnergy = solnData(ener_map,i,j,k) - kineticEnergy
                 end if
              end if
           else if(game_map /= NONEXISTENT) then ! This case should be usable for R(elativistic)HD - KW
              internalEnergy  = solnData(pres_map,i,j,k) / solnData(dens_map,i,j,k) / &
                                   (solnData(game_map,i,j,k) - 1.0)
              if(ener_map /= NONEXISTENT) then
                 if ( solnData(ener_map,i,j,k) - kineticEnergy > max(eos_smalle, eos_eintSwitch*kineticEnergy)) then
                    internalEnergy = solnData(ener_map,i,j,k) - kineticEnergy
                 end if
              end if
           else if(ener_map /= NONEXISTENT) then
              internalEnergy = solnData(ener_map,i,j,k)-kineticEnergy
           else
              internalEnergy = eos_smalle
           endif
           
           internalEnergy = max(internalEnergy, eos_smalle)
           eosData(eint+n) = internalEnergy
           
           eosData(pres+n) = solnData(pres_map,i,j,k)
           eosData(dens+n) = solnData(dens_map,i,j,k)
           eosData(temp+n) = solnData(temp_map,i,j,k)
           eosData(gamc+n) = solnData(gamc_map,i,j,k)
           if((ye_map /= NONEXISTENT).and.(sumy_map /= NONEXISTENT)) then
              !! cal says abar=1/sumy
              !! cal says zbar=ye / sumy and he claims sumy are never zero
              eosData(abar+n) =  1.0 /  solnData(sumy_map,i,j,k)
              eosData(zbar+n) = solnData(ye_map,i,j,k) /  solnData(sumy_map,i,j,k)
           endif
           if(entr_map /= NONEXISTENT) eosData(entr+n) = solnData(entr_map,i,j,k)
           
end subroutine Eos_getData_loop1