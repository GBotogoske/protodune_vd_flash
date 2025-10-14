#!/bin/bash

hv_list=("39516" "39517" "39518" "39519" "39521" "39522" "39523" "39525" "39526" "39528" "39529")

for i in "${hv_list[@]}"
do
    echo ">>> Processando run $i"
    
    ../build/flash -run "$i"
    ../build/filter_flash -run "$i"
    
    # Remove arquivos *_filtered.root que contenham o número da run
    rm /home/gabriel/Documents/protodune/data/VD/*"${i}"_flash.root
    
    ../build//avg -run "$i"
    
    echo ">>> Run $i concluída"
    echo "-----------------------------"
done
