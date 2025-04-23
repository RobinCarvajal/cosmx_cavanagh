#!/bin/bash

de_ant_condition="/mnt/f/cosmx_scripts/regional/differential_expression/condition/de_ant_condition.sh"
de_pos_condition="/mnt/f/cosmx_scripts/regional/differential_expression/condition/de_pos_condition.sh"

bash "${de_ant_condition}"
bash "${de_pos_condition}"
