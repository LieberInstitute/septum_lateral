#!/bin/bash

MAINDIR="/dcs04/lieber/marmaypag/pilotLS_LIBD1070"

find ${MAINDIR} -type d -exec nfs4_setfacl -a "A:g:hickslab@cm.cluster:RWX" {} \;
find ${MAINDIR} -type d -exec nfs4_setfacl -a "A:gfdi:hickslab@cm.cluster:RWX" {} \;
find ${MAINDIR} -type f -exec nfs4_setfacl -a "A:g:hickslab@cm.cluster:RW" {} \;

find ${MAINDIR} -type d -exec nfs4_setfacl -a "A:g:lieber_lcolladotor@cm.cluster:RWX" {} \;
find ${MAINDIR} -type d -exec nfs4_setfacl -a "A:gfdi:lieber_lcolladotor@cm.cluster:RWX" {} \;
find ${MAINDIR} -type f -exec nfs4_setfacl -a "A:g:lieber_lcolladotor@cm.cluster:RW" {} \;

find ${MAINDIR} -type d -exec nfs4_setfacl -a "A:g:lieber_marmaypag@cm.cluster:RWX" {} \;
find ${MAINDIR} -type d -exec nfs4_setfacl -a "A:gfdi:lieber_marmaypag@cm.cluster:RWX" {} \;
find ${MAINDIR} -type f -exec nfs4_setfacl -a "A:g:lieber_marmaypag@cm.cluster:RW" {} \;

## To move away from lieber_jaffe
chgrp lieber_marmaypag -R ${MAINDIR}

## For setting the group sticky bit
find ${MAINDIR} -type d | xargs chmod g+s

## Check settings
nfs4_getfacl ${MAINDIR}
