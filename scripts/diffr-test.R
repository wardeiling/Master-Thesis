library(diffr)

# Load the two files to compare
diffr::diffr(file2 = "qian2020/updated scripts/generative_model_updated_newmodels.R", file1 = "qian2020/original scripts/generative_model.R", 
             minJumpSize = 100)  # wordWrap = 100)

diffr::diffr(file2 = "qian2020/updated scripts/simulation_updated_newmodels.R", file1 = "qian2020/original scripts/simulation.R", 
             minJumpSize = 100)  # wordWrap = 100)
