#!/usr/bin/env bash

# Symlink this script from the build folder

make -j 8 

TEST_FOLDER=platforms/reference/tests
eval $TEST_FOLDER/TestReferenceMBPolElectrostaticsForce
eval $TEST_FOLDER/TestReferenceMBPolOneBodyForce
eval $TEST_FOLDER/TestReferenceMBPolTwoBodyForce
eval $TEST_FOLDER/TestReferenceMBPolThreeBodyForce
eval $TEST_FOLDER/TestReferenceMBPolDispersionForce
eval $TEST_FOLDER/TestReferenceMBPolIntegrationTest
eval $TEST_FOLDER/TestReferenceMBPol14WaterTest
