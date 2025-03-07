function b=verifyresidual_bvp5c(maxres)

global OCMATCONT OCBVP

b=(OCBVP.solutionVerificationStage && maxres < OCMATCONT.OPTIONS.meshadaptreltol);