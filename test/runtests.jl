using Curvelet
using Base.Test

N           = 512
sample_data = randn(N,N)
plan        = Curvelet.planCurveletTransform(N, Curvelet.C1)
tmp         = Curvelet.curveletTransform(sample_data, plan)
sample_hat  = Curvelet.inverseCurveletTransform(tmp, plan)
@test norm(sample_hat.-sample_hat) < 0.01
