##### Create psmc class #####

psmc <- setClass("psmc", slots = list(name="character", truncated="character", theta="numeric", rho="numeric", gen="numeric",
                                      u="numeric", pattern="character", n_free_lambdas="numeric", n_iterations="numeric", 
                                      N0="numeric", lk0="numeric", base_times="numeric", base_nes="numeric", times="numeric", 
                                      nes="numeric", ms_theta="numeric", ms_rho="numeric", ms_times = "numeric", ms_nes = "numeric", 
                                      ms_cmd="character", macs_theta="numeric", macs_rho="numeric", macs_cmd="character"))


