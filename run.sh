#!/bin/bash

BASEDIR=dist
PREFERENCE_FILE=../input/example_pairs/file_of_preferences.txt
MIN_NB_USER_PER_RESOURCE=2
MAX_NB_USER_PER_RESOURCE=2
RESOURCE_DUPLICATE_MODE="FILE_BASED(../input/example_pairs/all_projects.csv,0,4)"
RESOURCE_OWNERSHIP_MODE="FILE_BASED(../input/example_pairs/all_projects.csv,0,3)"
PREFERENCE_MEANING="PERSONAL_INSATISFACTION"
OWNER_DESIRES="AT_LEAST_ONE_INSTANCE_PER_OWNER"
OUTPUT_MODE="LATEX_REPORT"

cd $BASEDIR

java -jar RessourceAllocator.jar "PREFERENCE_FILE:$PREFERENCE_FILE" "MIN_NB_USER_PER_RESOURCE:$MIN_NB_USER_PER_RESOURCE" "MAX_NB_USER_PER_RESOURCE:$MAX_NB_USER_PER_RESOURCE" "RESOURCE_DUPLICATE_MODE:$RESOURCE_DUPLICATE_MODE" "RESOURCE_OWNERSHIP_MODE:$RESOURCE_OWNERSHIP_MODE" "PREFERENCE_MEANING:$PREFERENCE_MEANING" "OWNER_DESIRES:$OWNER_DESIRES" "OUTPUT_MODE:$OUTPUT_MODE"