require("src/Curvelet.jl")
#using Curvelet
using Images
using Winston

N = 256
peppers = imread("circle.png")
sample_data = reshape(peppers[1,:,:],N,N)'*1.0
##
###sample_data = rand(1024,1024)
#win = Curvelet._gen_windows(size(sample_data,1), Curvelet.C1)
plan = Curvelet.planCurveletTransform(N, Curvelet.C1)
@time tmp = Curvelet.curveletTransform(sample_data, plan)
@time sample_hat  = Curvelet.inverseCurveletTransform(tmp, plan)
println("Done")
println("Error was: ", norm(sample_data-sample_hat))
println()

#
##Create an approximation
#
function vectorize_curvelet(cc::Matrix{Complex{Float64}})
    vec(cc)
end

function vectorize_curvelet(cc::Curvelet.CurveletCoeff)
    data = vectorize_curvelet(cc.LP)
    for i = cc.HP
        data = [data,vec(i)]
    end
    data
end

function build_sparse(image, plan, sparsity)
    tmp = Curvelet.curveletTransform(sample_data, plan)
    data = vectorize_curvelet(tmp)
    ordered = reverse(sort(abs(data)))
    curvelet_threshold(tmp, ordered[sparsity])
    Curvelet.inverseCurveletTransform(tmp, plan)
end
#for i=tmp.HP
#    data = [data vec(i)]
#end
#
#data = vec(data)
#ordered = reverse(sort(abs(data)))

#println(ordered[10])
#println(ordered[100])
#println(ordered[1000])
#println(ordered[10000])
#
#threshold = ordered[30000]

function curvelet_threshold(cc::Matrix{Complex{Float64}}, threshold)
    #cc[abs(cc) .< threshold] = 0
end

function curvelet_threshold(cc::Curvelet.CurveletCoeff, threshold)
    curvelet_threshold(cc.LP, threshold)
    for i = 1:length(cc.HP)
        cc.HP[i][abs(cc.HP[i]) .< threshold] = 0
    end
end

#attempt a reconstruction problem
#grab around 50% of the points
observation = sample_data
for i=1:N, j=1:N
    observation[i,j] = (rand() < 0.2) ? sample_data[i,j] : 0
end

K=20000

#curvelet_threshold(tmp, threshold)
#@time sample_hat  = Curvelet.inverseCurveletTransform(tmp, plan)
#println()
#println()
#println()
#println("Appx error was: ", norm(sample_data-sample_hat))
#
#imagesc(abs(sample_hat))
#
#sleep(60)
#
#for i=length(tmp.HP)
#    tmp.HP[i][abs(tmp.HP[i]) .< ordered[10000]] = 0
#end
#
##tmp.LP[abs(tmp.LP) .< ordered[10000]] = 0
#
#@time sample_hat  = inverseCurveletTransform(tmp, plan)
#println("total energy:   ", norm(sample_data))
#println("Appx error was: ", norm(sample_data-sample_hat))
