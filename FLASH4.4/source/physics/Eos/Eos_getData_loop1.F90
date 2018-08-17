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

return          
end subroutine Eos_getData_loop1
