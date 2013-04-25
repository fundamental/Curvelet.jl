require("src/Curvelet.jl")
using Curvelet
using Images

N = 512
peppers = imread("test/peppers.gif")
sample_data = reshape(peppers[1,:,:],N,N)'*1.0
##
###sample_data = rand(1024,1024)
win = Curvelet._gen_windows(size(sample_data,1), Curvelet.C1)
plan = Curvelet.planCurveletTransform(N, Curvelet.C1)
@time tmp = curveletTransform(sample_data, plan)
@time sample_hat  = inverseCurveletTransform(tmp, plan)
println("Done")
println("Error was: ", norm(sample_data-sample_hat))
println()
