include("../src/Curvelet.jl")
#using Curvelet
using PyPlot

#Load sample data
input = imread("test/peppers.gif")
input = imread("test/lena.jpg")
#input = imread("smile.png")
input = input[:,:,1]
sample_data = input*1.0
N = size(sample_data,1)

#Display unmodified data
figure(1)
gray()
imshow(sample_data)
title("Unmodified Data")


plan = Curvelet.planCurveletTransform(N, Curvelet.C1)
@time tmp = Curvelet.curveletTransform(sample_data, plan)
@time sample_hat  = Curvelet.inverseCurveletTransform(tmp, plan)
println("Done")
println("Error from transform was: ", norm(sample_data-sample_hat))
println()
combo = zeros(size(tmp.HP[1]))
N = length(tmp.HP)
for i=1:2:length(tmp.HP)
    figure(i)
    t    = norm.(tmp.HP[i])
    combo += t
    imshow(log.(t))
end
figure(N+1)
imshow(log.(combo))

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
data = Float64[0.0]
for i=tmp.HP
    data = vcat(data, vec(i))
end
#
data = vec(data)
ordered = reverse(sort(abs(data)))

#threshold = ordered[30000]
threshold = ordered[300]

function curvelet_threshold(cc::Matrix{Complex{Float64}}, threshold)
    #cc[abs(cc) .< threshold] = 0
end

function curvelet_threshold(cc::Curvelet.CurveletCoeff, threshold)
    curvelet_threshold(cc.LP, threshold)
    for i = 1:length(cc.HP)
        cc.HP[i][abs.(cc.HP[i]) .< threshold] = 0
    end
end

#attempt a reconstruction problem
#grab around 75% of the points
observation = copy(sample_data)
N = size(observation,1)
for i=1:N, j=1:N
    observation[i,j] = (rand() < 0.75) ? sample_data[i,j] : 0
end

@time tmp = Curvelet.curveletTransform(observation, plan)
K=20000
data = Float64[0.0]
for i=tmp.HP
    data = vcat(data, vec(i))
end
#
data = vec(data)
ordered = reverse(sort(abs(data)))
threshold = ordered[3000]

curvelet_threshold(tmp, threshold)
@time sample_hat  = Curvelet.inverseCurveletTransform(tmp, plan)
println()
println()
println()
println("Appx error  was: ", norm(sample_data-sample_hat))
println("Orig energy was: ", norm(sample_data))

figure(100)
title("Approximation Under Missing Data - Src Data")
gray()
imshow(abs.(sample_data))
figure(101)
title("Approximation Under Missing Data - Missing")
gray()
imshow(abs.(observation))
figure(102)
title("Approximation Under Missing Data - Reconstruct")
gray()
imshow(abs.(sample_hat))


@time tmp = Curvelet.curveletTransform(sample_data, plan)
K = round(Int,length(sample_data)*0.10)
threshold = ordered[K]
curvelet_threshold(tmp,threshold)


@time sample_hat  = Curvelet.inverseCurveletTransform(tmp, plan)
println("Sparsity:       ", 100*K/length(sample_data), "%")
println("total energy:   ", norm(sample_data))
println("Appx error was: ", norm(sample_data-sample_hat))

figure(200)
title("Approximation Under Sparsity - Original")
gray()
imshow(abs.(sample_data))
figure(201)
title(string("Approximation Under Sparsity - ",K," sparse"))
gray()
imshow(abs.(sample_hat))
