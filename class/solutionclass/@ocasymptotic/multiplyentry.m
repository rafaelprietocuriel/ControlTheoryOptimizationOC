function ocAsym=multiplyentry(ocAsym,multiplyposition,multiple)

ocAsym=ocasymptotic(multiplyentry(octrajectory(ocAsym),multiplyposition,multiple),ocAsym.dynprimitive);